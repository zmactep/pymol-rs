//! Ring-buffer history of recent frame times for HUD percentile display.
//!
//! `FrameStats` from `crate::stats` is a single fully-resolved frame.
//! The HUD needs a *window* of recent frames so it can surface average
//! fps plus tail latency (1 % / 0.1 % low) — averages alone hide
//! micro-stutters.
//!
//! This module is compiled unconditionally — even when the `stats`
//! feature is off, the host can push wall-clock frame-gap times here to
//! get fps + percentile data without GPU timing.

use std::collections::VecDeque;

/// Default window size. Plan calls for 120-240 frames; 240 ≈ 2 s at
/// 120 fps which is enough to surface stutters without lag in the HUD.
pub const DEFAULT_CAPACITY: usize = 240;

/// Ring buffer of recent frame durations in milliseconds.
#[derive(Debug, Clone)]
pub struct FrameStatsHistory {
    ring: VecDeque<f32>,
    capacity: usize,
}

impl FrameStatsHistory {
    pub fn new(capacity: usize) -> Self {
        let cap = capacity.max(1);
        Self {
            ring: VecDeque::with_capacity(cap),
            capacity: cap,
        }
    }

    pub fn with_default_capacity() -> Self {
        Self::new(DEFAULT_CAPACITY)
    }

    pub fn push(&mut self, frame_ms: f32) {
        if frame_ms < 0.0 {
            return;
        }
        if self.ring.len() == self.capacity {
            self.ring.pop_front();
        }
        self.ring.push_back(frame_ms);
    }

    pub fn len(&self) -> usize {
        self.ring.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ring.is_empty()
    }

    pub fn clear(&mut self) {
        self.ring.clear();
    }

    /// Arithmetic mean over the window. Returns `0.0` on empty.
    pub fn avg_ms(&self) -> f32 {
        if self.ring.is_empty() {
            return 0.0;
        }
        let sum: f32 = self.ring.iter().sum();
        sum / self.ring.len() as f32
    }

    /// Average fps derived from the average frame time.
    pub fn avg_fps(&self) -> f32 {
        let m = self.avg_ms();
        if m <= 0.0 {
            0.0
        } else {
            1000.0 / m
        }
    }

    /// `p`-th percentile frame time, `p` in `[0.0, 1.0]`. Returns `0.0`
    /// on empty. O(n log n) per call (sorted copy of the window); the
    /// window is at most a few hundred entries so this is fine for HUD
    /// refresh frequency.
    pub fn percentile_ms(&self, p: f32) -> f32 {
        if self.ring.is_empty() {
            return 0.0;
        }
        let mut tmp: Vec<f32> = self.ring.iter().copied().collect();
        tmp.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let p = p.clamp(0.0, 1.0);
        let idx = ((tmp.len() as f32 - 1.0) * p).round() as usize;
        tmp[idx.min(tmp.len() - 1)]
    }

    /// "1 % low" — 95th-percentile frame time (= 5 %-of-frames-worse
    /// latency). This surfaces stutter alongside avg fps.
    pub fn low_1pct_ms(&self) -> f32 {
        self.percentile_ms(0.95)
    }

    /// "0.1 % low" — 99.9th-percentile frame time. Worst-case spikes.
    pub fn low_0_1pct_ms(&self) -> f32 {
        self.percentile_ms(0.999)
    }
}

impl Default for FrameStatsHistory {
    fn default() -> Self {
        Self::with_default_capacity()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_history_returns_zeros() {
        let h = FrameStatsHistory::new(10);
        assert_eq!(h.avg_ms(), 0.0);
        assert_eq!(h.avg_fps(), 0.0);
        assert_eq!(h.low_1pct_ms(), 0.0);
        assert_eq!(h.low_0_1pct_ms(), 0.0);
    }

    #[test]
    fn ring_evicts_oldest_when_full() {
        let mut h = FrameStatsHistory::new(3);
        h.push(1.0);
        h.push(2.0);
        h.push(3.0);
        assert_eq!(h.len(), 3);
        h.push(4.0); // evicts 1.0
        assert_eq!(h.len(), 3);
        // Sum = 2+3+4 = 9 → avg 3.0.
        assert!((h.avg_ms() - 3.0).abs() < 1e-5);
    }

    #[test]
    fn percentile_picks_extremes() {
        let mut h = FrameStatsHistory::new(100);
        for i in 1..=100 {
            h.push(i as f32);
        }
        // 95th percentile of 1..=100 ≈ 95.
        assert!((h.percentile_ms(0.95) - 95.0).abs() < 1e-3);
        // 50th percentile ≈ median ≈ 50.
        assert!((h.percentile_ms(0.5) - 51.0).abs() < 2.0);
        // 0th percentile = min.
        assert!((h.percentile_ms(0.0) - 1.0).abs() < 1e-3);
        // 100th percentile = max.
        assert!((h.percentile_ms(1.0) - 100.0).abs() < 1e-3);
    }

    #[test]
    fn negative_frames_ignored() {
        let mut h = FrameStatsHistory::new(10);
        h.push(-1.0);
        h.push(5.0);
        h.push(-2.0);
        assert_eq!(h.len(), 1);
        assert!((h.avg_ms() - 5.0).abs() < 1e-5);
    }

    #[test]
    fn avg_fps_from_frame_time() {
        let mut h = FrameStatsHistory::new(10);
        for _ in 0..10 {
            h.push(8.333); // ≈ 120 fps
        }
        let fps = h.avg_fps();
        assert!((fps - 120.0).abs() < 0.1, "got {fps}");
    }
}
