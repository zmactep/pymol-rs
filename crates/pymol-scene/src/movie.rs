//! Movie system for frame-based animation
//!
//! Provides keyframe animation and playback control for scenes.
//! Supports view interpolation, object state changes, and command execution.

use ahash::AHashMap;
use std::time::{Duration, Instant};

use crate::camera::SceneView;

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
        self.current_frame = 0;
        self.stop();
    }

    /// Set the number of frames (resizing the movie)
    pub fn set_frame_count(&mut self, count: usize) {
        self.frames.resize_with(count, MovieFrame::new);
        if self.current_frame >= count && count > 0 {
            self.current_frame = count - 1;
        }
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

    /// Get the interpolated view between current and next frame
    pub fn interpolated_view(&self) -> Option<SceneView> {
        if !self.interpolate || self.frames.is_empty() {
            return self.frames.get(self.current_frame).and_then(|f| f.view.clone());
        }

        let current_view = self.frames.get(self.current_frame).and_then(|f| f.view.clone());
        
        // Find next frame with a view
        let next_idx = match self.direction {
            PlayDirection::Forward => {
                if self.current_frame + 1 < self.frames.len() {
                    self.current_frame + 1
                } else {
                    return current_view;
                }
            }
            PlayDirection::Backward => {
                if self.current_frame > 0 {
                    self.current_frame - 1
                } else {
                    return current_view;
                }
            }
        };

        let next_view = self.frames.get(next_idx).and_then(|f| f.view.clone());

        match (current_view, next_view) {
            (Some(current), Some(next)) => {
                Some(Self::lerp_view(&current, &next, self.interpolation_t))
            }
            (Some(current), None) => Some(current),
            (None, Some(next)) => Some(next),
            (None, None) => None,
        }
    }

    /// Linear interpolation between two views
    fn lerp_view(from: &SceneView, to: &SceneView, t: f32) -> SceneView {
        let t = t.clamp(0.0, 1.0);
        
        // Interpolate position
        let position = lin_alg::f32::Vec3::new(
            from.position.x + (to.position.x - from.position.x) * t,
            from.position.y + (to.position.y - from.position.y) * t,
            from.position.z + (to.position.z - from.position.z) * t,
        );

        // Interpolate origin
        let origin = lin_alg::f32::Vec3::new(
            from.origin.x + (to.origin.x - from.origin.x) * t,
            from.origin.y + (to.origin.y - from.origin.y) * t,
            from.origin.z + (to.origin.z - from.origin.z) * t,
        );

        // Interpolate clip planes and FOV
        let clip_front = from.clip_front + (to.clip_front - from.clip_front) * t;
        let clip_back = from.clip_back + (to.clip_back - from.clip_back) * t;
        let fov = from.fov + (to.fov - from.fov) * t;

        // For rotation, we should use quaternion slerp, but for simplicity
        // we'll interpolate the rotation matrix directly (works for small angles)
        // A proper implementation would extract quaternions and slerp them
        let rotation = Self::lerp_mat4(&from.rotation, &to.rotation, t);

        SceneView {
            rotation,
            position,
            origin,
            clip_front,
            clip_back,
            fov,
        }
    }

    /// Linear interpolation between matrices (not ideal for rotations, but works for small changes)
    fn lerp_mat4(from: &lin_alg::f32::Mat4, to: &lin_alg::f32::Mat4, t: f32) -> lin_alg::f32::Mat4 {
        let mut result = lin_alg::f32::Mat4::new_identity();
        for i in 0..16 {
            result.data[i] = from.data[i] + (to.data[i] - from.data[i]) * t;
        }
        result
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
