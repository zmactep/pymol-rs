//! Movie system for frame-based animation
//!
//! Provides keyframe animation and playback control for scenes.
//! Supports view interpolation, object state changes, and command execution.

use ahash::AHashMap;
use std::time::{Duration, Instant};

use crate::camera::SceneView;
use crate::quat::{self, Quat};

/// Movie loop mode
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LoopMode {
    /// Play once and stop
    #[default]
    Once,
    /// Loop continuously
    Loop,
    /// Play forward then backward (ping-pong)
    Swing,
}

/// Playback direction
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlayDirection {
    Forward,
    Backward,
}

impl Default for PlayDirection {
    fn default() -> Self {
        Self::Forward
    }
}

/// Object-specific keyframe data
#[derive(Debug, Clone)]
pub struct ObjectKeyframe {
    /// Object transform matrix
    pub transform: Option<lin_alg::f32::Mat4>,
    /// Coordinate set state
    pub state: Option<usize>,
}

/// A single movie frame
#[derive(Debug, Clone, Default)]
pub struct MovieFrame {
    /// Camera view for this frame (if stored)
    pub view: Option<SceneView>,
    /// Object states: object name -> coordinate state index
    pub object_states: AHashMap<String, usize>,
    /// Commands to execute when reaching this frame
    pub commands: Vec<String>,
    /// Message to display
    pub message: Option<String>,
    /// Scene keyframe reference (name of a stored scene)
    pub scene_name: Option<String>,
    /// Whether this frame is a camera keyframe
    pub is_camera_keyframe: bool,
    /// Per-object keyframes
    pub object_keyframes: AHashMap<String, ObjectKeyframe>,
    /// Global state index for this frame (from mset specification)
    pub state_index: usize,
}

impl MovieFrame {
    /// Create a new empty frame
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a frame with a view
    pub fn with_view(view: SceneView) -> Self {
        Self {
            view: Some(view),
            ..Default::default()
        }
    }

    /// Set the view for this frame
    pub fn set_view(&mut self, view: SceneView) {
        self.view = Some(view);
    }

    /// Set an object's state for this frame
    pub fn set_object_state(&mut self, object_name: &str, state: usize) {
        self.object_states.insert(object_name.to_string(), state);
    }

    /// Add a command to execute
    pub fn add_command(&mut self, command: &str) {
        self.commands.push(command.to_string());
    }

    /// Check if this frame has any data
    pub fn is_empty(&self) -> bool {
        self.view.is_none() && self.object_states.is_empty() && self.commands.is_empty()
    }
}

/// Movie playback state
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlaybackState {
    /// Movie is stopped at a specific frame
    Stopped,
    /// Movie is playing
    Playing,
    /// Movie is paused
    Paused,
}

impl Default for PlaybackState {
    fn default() -> Self {
        Self::Stopped
    }
}

/// Movie player for frame-based animation
///
/// Manages a sequence of frames with views, object states, and commands.
/// Supports various playback modes including looping and swing (ping-pong).
#[derive(Debug)]
pub struct Movie {
    /// Movie frames
    frames: Vec<MovieFrame>,
    /// Current frame index
    current_frame: usize,
    /// Playback state
    state: PlaybackState,
    /// Loop mode
    loop_mode: LoopMode,
    /// Playback direction (for swing mode)
    direction: PlayDirection,
    /// Frames per second
    fps: f32,
    /// Frame delay (derived from fps)
    frame_delay: Duration,
    /// Time of last frame advance
    last_frame_time: Option<Instant>,
    /// Whether to interpolate views between keyframes
    interpolate: bool,
    /// Interpolation progress within current frame (0.0 to 1.0)
    interpolation_t: f32,
    /// Cache for rock/wobble motion
    rock_angle: f32,
    /// Whether rock mode is enabled
    rock_enabled: bool,
    /// Precomputed interpolated views for all frames between keyframes
    precomputed_views: Vec<Option<SceneView>>,
}

impl Default for Movie {
    fn default() -> Self {
        Self::new()
    }
}

impl Movie {
    /// Create a new empty movie
    pub fn new() -> Self {
        Self {
            frames: Vec::new(),
            current_frame: 0,
            state: PlaybackState::Stopped,
            loop_mode: LoopMode::default(),
            direction: PlayDirection::Forward,
            fps: 30.0,
            frame_delay: Duration::from_secs_f32(1.0 / 30.0),
            last_frame_time: None,
            interpolate: true,
            interpolation_t: 0.0,
            rock_angle: 0.0,
            rock_enabled: false,
            precomputed_views: Vec::new(),
        }
    }

    /// Create a movie with a specific number of empty frames
    pub fn with_frames(count: usize) -> Self {
        let mut movie = Self::new();
        movie.frames = vec![MovieFrame::new(); count];
        movie
    }

    /// Get the number of frames
    pub fn frame_count(&self) -> usize {
        self.frames.len()
    }

    /// Check if the movie is empty
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// Get the current frame index
    pub fn current_frame(&self) -> usize {
        self.current_frame
    }

    /// Get the current frame
    pub fn frame(&self) -> Option<&MovieFrame> {
        self.frames.get(self.current_frame)
    }

    /// Get a frame by index
    pub fn frame_at(&self, index: usize) -> Option<&MovieFrame> {
        self.frames.get(index)
    }

    /// Get mutable access to a frame
    pub fn frame_mut(&mut self, index: usize) -> Option<&mut MovieFrame> {
        self.frames.get_mut(index)
    }

    /// Add a frame to the end of the movie
    pub fn add_frame(&mut self, frame: MovieFrame) {
        self.frames.push(frame);
    }

    /// Insert a frame at a specific index
    pub fn insert_frame(&mut self, index: usize, frame: MovieFrame) {
        if index <= self.frames.len() {
            self.frames.insert(index, frame);
        }
    }

    /// Remove a frame
    pub fn remove_frame(&mut self, index: usize) -> Option<MovieFrame> {
        if index < self.frames.len() {
            let frame = self.frames.remove(index);
            // Adjust current frame if necessary
            if self.current_frame >= self.frames.len() && !self.frames.is_empty() {
                self.current_frame = self.frames.len() - 1;
            }
            Some(frame)
        } else {
            None
        }
    }

    /// Clear all frames
    pub fn clear(&mut self) {
        self.frames.clear();
        self.precomputed_views.clear();
        self.current_frame = 0;
        self.stop();
    }

    /// Set the number of frames (resizing the movie)
    pub fn set_frame_count(&mut self, count: usize) {
        self.frames.resize_with(count, MovieFrame::new);
        self.precomputed_views.clear();
        if self.current_frame >= count && count > 0 {
            self.current_frame = count - 1;
        }
    }

    /// Set frames from an mset specification (state index per frame).
    ///
    /// The `states` vec contains one state index per frame.
    pub fn set_from_spec(&mut self, states: Vec<usize>) {
        let count = states.len();
        self.frames.resize_with(count, MovieFrame::new);
        for (i, state) in states.into_iter().enumerate() {
            self.frames[i].state_index = state;
        }
        self.precomputed_views.clear();
        if self.current_frame >= count && count > 0 {
            self.current_frame = count - 1;
        }
    }

    /// Append frames from an mset specification onto the existing movie.
    pub fn append_from_spec(&mut self, states: Vec<usize>) {
        for state in states {
            let mut frame = MovieFrame::new();
            frame.state_index = state;
            self.frames.push(frame);
        }
        self.precomputed_views.clear();
    }

    /// Get the FPS
    pub fn fps(&self) -> f32 {
        self.fps
    }

    /// Set the FPS
    pub fn set_fps(&mut self, fps: f32) {
        self.fps = fps.max(0.1);
        self.frame_delay = Duration::from_secs_f32(1.0 / self.fps);
    }

    /// Get the loop mode
    pub fn loop_mode(&self) -> LoopMode {
        self.loop_mode
    }

    /// Set the loop mode
    pub fn set_loop_mode(&mut self, mode: LoopMode) {
        self.loop_mode = mode;
    }

    /// Check if interpolation is enabled
    pub fn interpolate(&self) -> bool {
        self.interpolate
    }

    /// Set whether to interpolate between frames
    pub fn set_interpolate(&mut self, interpolate: bool) {
        self.interpolate = interpolate;
    }

    /// Get the playback state
    pub fn playback_state(&self) -> PlaybackState {
        self.state
    }

    /// Check if playing
    pub fn is_playing(&self) -> bool {
        self.state == PlaybackState::Playing
    }

    /// Check if paused
    pub fn is_paused(&self) -> bool {
        self.state == PlaybackState::Paused
    }

    /// Check if stopped
    pub fn is_stopped(&self) -> bool {
        self.state == PlaybackState::Stopped
    }

    /// Start playing
    pub fn play(&mut self) {
        if !self.frames.is_empty() {
            self.state = PlaybackState::Playing;
            self.last_frame_time = Some(Instant::now());
        }
    }

    /// Pause playback
    pub fn pause(&mut self) {
        if self.state == PlaybackState::Playing {
            self.state = PlaybackState::Paused;
        }
    }

    /// Resume from pause
    pub fn resume(&mut self) {
        if self.state == PlaybackState::Paused {
            self.state = PlaybackState::Playing;
            self.last_frame_time = Some(Instant::now());
        }
    }

    /// Stop playback and reset to first frame
    pub fn stop(&mut self) {
        self.state = PlaybackState::Stopped;
        self.current_frame = 0;
        self.direction = PlayDirection::Forward;
        self.interpolation_t = 0.0;
        self.last_frame_time = None;
    }

    /// Toggle play/pause
    pub fn toggle(&mut self) {
        match self.state {
            PlaybackState::Playing => self.pause(),
            PlaybackState::Paused => self.resume(),
            PlaybackState::Stopped => self.play(),
        }
    }

    /// Go to a specific frame
    pub fn goto_frame(&mut self, frame: usize) {
        if frame < self.frames.len() {
            self.current_frame = frame;
            self.interpolation_t = 0.0;
        }
    }

    /// Go to the first frame
    pub fn rewind(&mut self) {
        self.goto_frame(0);
        self.direction = PlayDirection::Forward;
    }

    /// Go to the last frame
    pub fn fast_forward(&mut self) {
        if !self.frames.is_empty() {
            self.goto_frame(self.frames.len() - 1);
        }
    }

    /// Advance to next frame (manual)
    pub fn next_frame(&mut self) {
        if self.current_frame + 1 < self.frames.len() {
            self.current_frame += 1;
            self.interpolation_t = 0.0;
        }
    }

    /// Go to previous frame (manual)
    pub fn prev_frame(&mut self) {
        if self.current_frame > 0 {
            self.current_frame -= 1;
            self.interpolation_t = 0.0;
        }
    }

    /// Update the movie, advancing frames based on time
    ///
    /// Returns true if the frame changed
    pub fn update(&mut self) -> bool {
        if self.state != PlaybackState::Playing || self.frames.is_empty() {
            return false;
        }

        let now = Instant::now();
        let elapsed = self.last_frame_time.map_or(Duration::ZERO, |t| now - t);

        if elapsed < self.frame_delay {
            // Update interpolation progress
            if self.interpolate {
                self.interpolation_t = elapsed.as_secs_f32() / self.frame_delay.as_secs_f32();
            }
            return false;
        }

        // Time to advance frame
        self.last_frame_time = Some(now);
        self.interpolation_t = 0.0;

        let frame_changed = self.advance_frame();
        frame_changed
    }

    /// Advance to the next frame based on direction and loop mode
    fn advance_frame(&mut self) -> bool {
        let old_frame = self.current_frame;

        match self.direction {
            PlayDirection::Forward => {
                if self.current_frame + 1 < self.frames.len() {
                    self.current_frame += 1;
                } else {
                    // Reached end
                    match self.loop_mode {
                        LoopMode::Once => {
                            self.state = PlaybackState::Stopped;
                        }
                        LoopMode::Loop => {
                            self.current_frame = 0;
                        }
                        LoopMode::Swing => {
                            self.direction = PlayDirection::Backward;
                            if self.current_frame > 0 {
                                self.current_frame -= 1;
                            }
                        }
                    }
                }
            }
            PlayDirection::Backward => {
                if self.current_frame > 0 {
                    self.current_frame -= 1;
                } else {
                    // Reached beginning
                    match self.loop_mode {
                        LoopMode::Once => {
                            self.state = PlaybackState::Stopped;
                        }
                        LoopMode::Loop => {
                            self.current_frame = self.frames.len().saturating_sub(1);
                        }
                        LoopMode::Swing => {
                            self.direction = PlayDirection::Forward;
                            if self.current_frame + 1 < self.frames.len() {
                                self.current_frame += 1;
                            }
                        }
                    }
                }
            }
        }

        self.current_frame != old_frame
    }

    /// Get the view for the current frame.
    ///
    /// Prefers precomputed interpolated views. Falls back to per-frame views
    /// or simple linear interpolation between adjacent frames.
    pub fn interpolated_view(&self) -> Option<SceneView> {
        if self.frames.is_empty() {
            return None;
        }

        // Use precomputed views if available
        if let Some(view) = self.precomputed_views.get(self.current_frame).and_then(|v| v.clone()) {
            return Some(view);
        }

        // Fall back to per-frame stored view
        self.frames.get(self.current_frame).and_then(|f| f.view.clone())
    }

    // =========================================================================
    // Keyframe Management
    // =========================================================================

    /// Store a camera keyframe at the given frame index.
    pub fn store_camera_keyframe(&mut self, frame: usize, view: SceneView) {
        if frame < self.frames.len() {
            self.frames[frame].view = Some(view);
            self.frames[frame].is_camera_keyframe = true;
        }
    }

    /// Store a scene keyframe at the given frame index.
    pub fn store_scene_keyframe(&mut self, frame: usize, scene_name: &str, view: SceneView) {
        if frame < self.frames.len() {
            self.frames[frame].scene_name = Some(scene_name.to_string());
            self.frames[frame].view = Some(view);
            self.frames[frame].is_camera_keyframe = true;
        }
    }

    /// Store an object keyframe at the given frame index.
    pub fn store_object_keyframe(&mut self, frame: usize, object: &str, state: Option<usize>) {
        if frame < self.frames.len() {
            self.frames[frame].object_keyframes.insert(
                object.to_string(),
                ObjectKeyframe {
                    transform: None,
                    state,
                },
            );
        }
    }

    /// Clear camera keyframe at a specific frame, or all if `frame` is None.
    pub fn clear_camera_keyframes(&mut self, frame: Option<usize>) {
        match frame {
            Some(f) => {
                if f < self.frames.len() {
                    self.frames[f].view = None;
                    self.frames[f].is_camera_keyframe = false;
                    self.frames[f].scene_name = None;
                }
            }
            None => {
                for f in &mut self.frames {
                    f.view = None;
                    f.is_camera_keyframe = false;
                    f.scene_name = None;
                }
            }
        }
        self.precomputed_views.clear();
    }

    /// Clear object keyframes at a specific frame for a specific object.
    pub fn clear_object_keyframes(&mut self, frame: Option<usize>, object: &str) {
        match frame {
            Some(f) => {
                if f < self.frames.len() {
                    self.frames[f].object_keyframes.remove(object);
                }
            }
            None => {
                for f in &mut self.frames {
                    f.object_keyframes.remove(object);
                }
            }
        }
    }

    // =========================================================================
    // Keyframe Interpolation
    // =========================================================================

    /// Compute smooth interpolated views for all frames between keyframes.
    ///
    /// Uses Catmull-Rom cubic Hermite splines for position/origin/clip/fov
    /// and quaternion SLERP for rotation.
    pub fn interpolate_keyframes(&mut self, loop_movie: bool) {
        let n = self.frames.len();
        if n == 0 {
            return;
        }

        // Collect camera keyframe indices and views
        let keyframes: Vec<(usize, SceneView)> = self
            .frames
            .iter()
            .enumerate()
            .filter(|(_, f)| f.is_camera_keyframe && f.view.is_some())
            .map(|(i, f)| (i, f.view.clone().unwrap()))
            .collect();

        if keyframes.is_empty() {
            self.precomputed_views = vec![None; n];
            return;
        }

        // Single keyframe: all frames get that view
        if keyframes.len() == 1 {
            let view = &keyframes[0].1;
            self.precomputed_views = vec![Some(view.clone()); n];
            return;
        }

        let mut views: Vec<Option<SceneView>> = vec![None; n];

        // Extract quaternions and scalar data for each keyframe
        let kf_quats: Vec<Quat> = keyframes
            .iter()
            .map(|(_, v)| Quat::from_mat4(&v.rotation).normalized())
            .collect();

        let kf_count = keyframes.len();

        // Interpolate between each consecutive pair of keyframes
        for seg in 0..kf_count {
            let next_seg = if seg + 1 < kf_count {
                seg + 1
            } else if loop_movie {
                0 // wrap around
            } else {
                break; // no more segments
            };

            let (start_frame, ref start_view) = keyframes[seg];
            let (end_frame, ref end_view) = keyframes[next_seg];

            // Compute span (number of frames in this segment)
            let span = if end_frame > start_frame {
                end_frame - start_frame
            } else if loop_movie && next_seg == 0 {
                // Wrap-around: frames from start to end of movie, then from 0 to end_frame
                n - start_frame + end_frame
            } else {
                continue;
            };

            if span == 0 {
                continue;
            }

            // Get tangent control points for Catmull-Rom
            let prev_seg = if seg > 0 {
                seg - 1
            } else if loop_movie {
                kf_count - 1
            } else {
                seg
            };
            let next_next_seg = if next_seg + 1 < kf_count {
                next_seg + 1
            } else if loop_movie {
                0
            } else {
                next_seg
            };

            let prev_view = &keyframes[prev_seg].1;
            let next_next_view = &keyframes[next_next_seg].1;

            // Catmull-Rom tangents for position
            let m0_pos = quat::catmull_rom_tangent_vec3(&prev_view.position, &end_view.position);
            let m1_pos =
                quat::catmull_rom_tangent_vec3(&start_view.position, &next_next_view.position);

            // Catmull-Rom tangents for origin
            let m0_org = quat::catmull_rom_tangent_vec3(&prev_view.origin, &end_view.origin);
            let m1_org =
                quat::catmull_rom_tangent_vec3(&start_view.origin, &next_next_view.origin);

            // Scalar tangents
            let m0_cf = quat::catmull_rom_tangent(prev_view.clip_front, end_view.clip_front);
            let m1_cf =
                quat::catmull_rom_tangent(start_view.clip_front, next_next_view.clip_front);
            let m0_cb = quat::catmull_rom_tangent(prev_view.clip_back, end_view.clip_back);
            let m1_cb = quat::catmull_rom_tangent(start_view.clip_back, next_next_view.clip_back);
            let m0_fov = quat::catmull_rom_tangent(prev_view.fov, end_view.fov);
            let m1_fov = quat::catmull_rom_tangent(start_view.fov, next_next_view.fov);

            let q0 = &kf_quats[seg];
            let q1 = &kf_quats[next_seg];

            for step in 0..=span {
                let frame_idx = (start_frame + step) % n;
                let t = step as f32 / span as f32;

                let position = quat::hermite_vec3(
                    &start_view.position,
                    &m0_pos,
                    &end_view.position,
                    &m1_pos,
                    t,
                );
                let origin = quat::hermite_vec3(
                    &start_view.origin,
                    &m0_org,
                    &end_view.origin,
                    &m1_org,
                    t,
                );
                let clip_front =
                    quat::hermite(start_view.clip_front, m0_cf, end_view.clip_front, m1_cf, t);
                let clip_back =
                    quat::hermite(start_view.clip_back, m0_cb, end_view.clip_back, m1_cb, t);
                let fov = quat::hermite(start_view.fov, m0_fov, end_view.fov, m1_fov, t);

                let rotation = Quat::slerp(q0, q1, t).to_mat4();

                views[frame_idx] = Some(SceneView {
                    rotation,
                    position,
                    origin,
                    clip_front,
                    clip_back,
                    fov,
                });
            }
        }

        // Fill frames before the first keyframe (hold first keyframe view)
        let first_kf_frame = keyframes[0].0;
        let first_view = keyframes[0].1.clone();
        for i in 0..first_kf_frame {
            if views[i].is_none() {
                views[i] = Some(first_view.clone());
            }
        }

        // Fill frames after the last keyframe (hold last keyframe view)
        if !loop_movie {
            let last_kf_frame = keyframes[kf_count - 1].0;
            let last_view = keyframes[kf_count - 1].1.clone();
            for i in (last_kf_frame + 1)..n {
                if views[i].is_none() {
                    views[i] = Some(last_view.clone());
                }
            }
        }

        self.precomputed_views = views;
    }

    /// Check if precomputed interpolation data is available.
    pub fn has_interpolation(&self) -> bool {
        !self.precomputed_views.is_empty()
    }

    /// Get the scene name for the current frame, if any.
    pub fn current_scene_name(&self) -> Option<&str> {
        self.frames
            .get(self.current_frame)
            .and_then(|f| f.scene_name.as_deref())
    }

    /// Enable or disable rock mode
    pub fn set_rock(&mut self, enabled: bool) {
        self.rock_enabled = enabled;
        self.rock_angle = 0.0;
    }

    /// Check if rock mode is enabled
    pub fn is_rock_enabled(&self) -> bool {
        self.rock_enabled
    }

    /// Get the current rock angle for view modification
    pub fn rock_angle(&self) -> f32 {
        self.rock_angle
    }

    /// Update rock animation
    ///
    /// Returns the delta angle to apply to the view
    pub fn update_rock(&mut self, dt: f32, amplitude: f32, speed: f32) -> f32 {
        if !self.rock_enabled {
            return 0.0;
        }

        self.rock_angle += dt * speed;
        amplitude * self.rock_angle.sin()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_movie_creation() {
        let movie = Movie::new();
        assert!(movie.is_empty());
        assert!(movie.is_stopped());
        assert_eq!(movie.fps(), 30.0);
    }

    #[test]
    fn test_movie_with_frames() {
        let movie = Movie::with_frames(10);
        assert_eq!(movie.frame_count(), 10);
        assert!(!movie.is_empty());
    }

    #[test]
    fn test_movie_frame_navigation() {
        let mut movie = Movie::with_frames(5);

        assert_eq!(movie.current_frame(), 0);
        movie.next_frame();
        assert_eq!(movie.current_frame(), 1);
        movie.prev_frame();
        assert_eq!(movie.current_frame(), 0);

        movie.goto_frame(3);
        assert_eq!(movie.current_frame(), 3);

        movie.fast_forward();
        assert_eq!(movie.current_frame(), 4);

        movie.rewind();
        assert_eq!(movie.current_frame(), 0);
    }

    #[test]
    fn test_movie_add_remove_frames() {
        let mut movie = Movie::new();

        movie.add_frame(MovieFrame::new());
        movie.add_frame(MovieFrame::new());
        assert_eq!(movie.frame_count(), 2);

        movie.remove_frame(0);
        assert_eq!(movie.frame_count(), 1);

        movie.clear();
        assert!(movie.is_empty());
    }

    #[test]
    fn test_movie_fps() {
        let mut movie = Movie::new();

        movie.set_fps(60.0);
        assert_eq!(movie.fps(), 60.0);

        movie.set_fps(0.0); // Should clamp to minimum
        assert!(movie.fps() >= 0.1);
    }

    #[test]
    fn test_movie_loop_mode() {
        let mut movie = Movie::new();

        assert_eq!(movie.loop_mode(), LoopMode::Once);
        movie.set_loop_mode(LoopMode::Loop);
        assert_eq!(movie.loop_mode(), LoopMode::Loop);
    }

    #[test]
    fn test_movie_playback_state() {
        let mut movie = Movie::with_frames(3);

        assert!(movie.is_stopped());

        movie.play();
        assert!(movie.is_playing());

        movie.pause();
        assert!(movie.is_paused());

        movie.resume();
        assert!(movie.is_playing());

        movie.stop();
        assert!(movie.is_stopped());
    }

    #[test]
    fn test_movie_toggle() {
        let mut movie = Movie::with_frames(3);

        movie.toggle(); // stopped -> playing
        assert!(movie.is_playing());

        movie.toggle(); // playing -> paused
        assert!(movie.is_paused());

        movie.toggle(); // paused -> playing
        assert!(movie.is_playing());
    }

    #[test]
    fn test_movie_frame_data() {
        let mut frame = MovieFrame::new();

        assert!(frame.is_empty());

        frame.set_object_state("protein", 5);
        frame.add_command("show cartoon");

        assert!(!frame.is_empty());
        assert_eq!(frame.object_states.get("protein"), Some(&5));
        assert_eq!(frame.commands.len(), 1);
    }

    #[test]
    fn test_movie_advance_loop() {
        let mut movie = Movie::with_frames(3);
        movie.set_loop_mode(LoopMode::Loop);
        movie.direction = PlayDirection::Forward;
        movie.current_frame = 2; // Last frame
        movie.state = PlaybackState::Playing;

        movie.advance_frame();
        assert_eq!(movie.current_frame(), 0); // Should wrap to first
    }

    #[test]
    fn test_movie_advance_swing() {
        let mut movie = Movie::with_frames(3);
        movie.set_loop_mode(LoopMode::Swing);
        movie.direction = PlayDirection::Forward;
        movie.current_frame = 2; // Last frame
        movie.state = PlaybackState::Playing;

        movie.advance_frame();
        assert_eq!(movie.current_frame(), 1); // Should reverse
        assert_eq!(movie.direction, PlayDirection::Backward);
    }
}
