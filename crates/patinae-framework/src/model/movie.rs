//! Movie timeline model.
//!
//! Converts the UI-agnostic [`patinae_scene::Movie`] state into compact lanes
//! that frontends can render without peeking into playback details.

use patinae_scene::{Movie, Session};

/// A ruler tick on the movie timeline.
#[derive(Debug, Clone, PartialEq)]
pub struct MovieTick {
    /// 1-based frame number.
    pub frame: usize,
    /// Display label.
    pub label: String,
    /// Center position along the timeline track, in the inclusive range 0..=1.
    pub position: f32,
}

/// A point marker on a timeline lane.
#[derive(Debug, Clone, PartialEq)]
pub struct MovieMarker {
    /// 1-based frame number.
    pub frame: usize,
    /// Center position along the timeline track, in the inclusive range 0..=1.
    pub position: f32,
    /// Marker kind: "camera", "scene", or "object".
    pub kind: String,
    /// Human-readable label.
    pub label: String,
}

/// A continuous timeline block.
#[derive(Debug, Clone, PartialEq)]
pub struct MovieBlock {
    /// 1-based first frame in the block.
    pub start_frame: usize,
    /// 1-based last frame in the block.
    pub end_frame: usize,
    /// Start cell edge, in the inclusive range 0..=1.
    pub start_position: f32,
    /// End cell edge, in the inclusive range 0..=1.
    pub end_position: f32,
    /// Block kind: "state" or "object-state".
    pub kind: String,
    /// Human-readable label.
    pub label: String,
}

/// A named timeline lane.
#[derive(Debug, Clone, PartialEq)]
pub struct MovieLane {
    /// Stable lane id.
    pub id: String,
    /// Short display label.
    pub label: String,
    /// Point markers on this lane.
    pub markers: Vec<MovieMarker>,
    /// Continuous blocks on this lane.
    pub blocks: Vec<MovieBlock>,
}

/// Complete compact movie timeline.
#[derive(Debug, Clone, PartialEq)]
pub struct MovieTimeline {
    /// Effective frame count: explicit movie frames or object-state count.
    pub frame_count: usize,
    /// Ruler ticks.
    pub ticks: Vec<MovieTick>,
    /// Renderable lanes.
    pub lanes: Vec<MovieLane>,
}

impl Default for MovieTimeline {
    fn default() -> Self {
        Self::from_movie(&Movie::new())
    }
}

impl MovieTimeline {
    /// Build a timeline from a movie.
    pub fn from_movie(movie: &Movie) -> Self {
        let frame_count = movie.effective_frame_count().max(1);
        let ticks = build_ticks(frame_count);

        let lanes = vec![
            build_camera_lane(movie, frame_count),
            build_scene_lane(movie, frame_count),
            build_object_lane(movie, frame_count),
            build_mset_lane(movie, frame_count),
        ];

        Self {
            frame_count,
            ticks,
            lanes,
        }
    }

    /// Find a lane by id.
    pub fn lane(&self, id: &str) -> Option<&MovieLane> {
        self.lanes.iter().find(|lane| lane.id == id)
    }
}

/// Cached timeline model for frontends that sync every frame.
#[derive(Debug, Default)]
pub struct MovieTimelineModel {
    timeline: MovieTimeline,
    last_movie_generation: Option<u64>,
    last_registry_generation: Option<u64>,
    last_frame_count: usize,
}

impl MovieTimelineModel {
    /// Current cached timeline.
    pub fn timeline(&self) -> &MovieTimeline {
        &self.timeline
    }

    /// Rebuild if the movie structure or object-state count changed.
    pub fn sync(&mut self, session: &Session) -> bool {
        let movie_generation = session.movie.generation();
        let registry_generation = session.registry.generation();
        let frame_count = session.movie.effective_frame_count();

        let changed = self.last_movie_generation != Some(movie_generation)
            || self.last_registry_generation != Some(registry_generation)
            || self.last_frame_count != frame_count;

        if changed {
            self.timeline = MovieTimeline::from_movie(&session.movie);
            self.last_movie_generation = Some(movie_generation);
            self.last_registry_generation = Some(registry_generation);
            self.last_frame_count = frame_count;
        }

        changed
    }
}

fn build_camera_lane(movie: &Movie, frame_count: usize) -> MovieLane {
    let mut markers = Vec::new();
    for i in 0..movie.frame_count() {
        let Some(frame) = movie.frame_at(i) else {
            continue;
        };
        if frame.is_camera_keyframe && frame.view.is_some() {
            markers.push(MovieMarker {
                frame: i + 1,
                position: frame_position(i, frame_count),
                kind: "camera".to_string(),
                label: "camera".to_string(),
            });
        }
    }

    MovieLane {
        id: "camera".to_string(),
        label: "CAM".to_string(),
        markers,
        blocks: Vec::new(),
    }
}

fn build_scene_lane(movie: &Movie, frame_count: usize) -> MovieLane {
    let mut markers = Vec::new();
    for i in 0..movie.frame_count() {
        let Some(frame) = movie.frame_at(i) else {
            continue;
        };
        if let Some(scene_name) = &frame.scene_name {
            markers.push(MovieMarker {
                frame: i + 1,
                position: frame_position(i, frame_count),
                kind: "scene".to_string(),
                label: scene_name.clone(),
            });
        }
    }

    MovieLane {
        id: "scene".to_string(),
        label: "SCN".to_string(),
        markers,
        blocks: Vec::new(),
    }
}

fn build_object_lane(movie: &Movie, frame_count: usize) -> MovieLane {
    let mut markers = Vec::new();
    let mut per_frame_labels = Vec::with_capacity(frame_count);

    for i in 0..frame_count {
        let Some(frame) = movie.frame_at(i) else {
            per_frame_labels.push(None);
            continue;
        };

        if !frame.object_keyframes.is_empty() {
            let label = if frame.object_keyframes.len() == 1 {
                frame
                    .object_keyframes
                    .keys()
                    .next()
                    .cloned()
                    .unwrap_or_else(|| "object".to_string())
            } else {
                format!("{} objects", frame.object_keyframes.len())
            };
            markers.push(MovieMarker {
                frame: i + 1,
                position: frame_position(i, frame_count),
                kind: "object".to_string(),
                label,
            });
        }

        let mut labels: Vec<String> = frame
            .object_states
            .iter()
            .map(|(name, state)| format!("{name}: state {state}"))
            .collect();
        labels.extend(
            frame
                .object_keyframes
                .iter()
                .filter_map(|(name, keyframe)| {
                    keyframe.state.map(|state| format!("{name}: state {state}"))
                }),
        );
        labels.sort();

        per_frame_labels.push(if labels.is_empty() {
            None
        } else {
            Some(labels.join(", "))
        });
    }

    MovieLane {
        id: "object".to_string(),
        label: "OBJ".to_string(),
        markers,
        blocks: build_optional_blocks(&per_frame_labels, frame_count, "object-state"),
    }
}

fn build_mset_lane(movie: &Movie, frame_count: usize) -> MovieLane {
    let labels: Vec<Option<String>> = (0..frame_count)
        .map(|frame| Some(format!("state {}", movie.frame_to_state(frame) + 1)))
        .collect();

    MovieLane {
        id: "mset".to_string(),
        label: "MSET".to_string(),
        markers: Vec::new(),
        blocks: build_optional_blocks(&labels, frame_count, "state"),
    }
}

fn build_optional_blocks(
    labels: &[Option<String>],
    frame_count: usize,
    kind: &str,
) -> Vec<MovieBlock> {
    let mut blocks = Vec::new();
    let mut start = 0usize;
    let mut current: Option<&str> = None;

    for i in 0..=labels.len() {
        let next = labels.get(i).and_then(|label| label.as_deref());
        if i == 0 {
            start = 0;
            current = next;
            continue;
        }

        if next != current {
            if let Some(label) = current {
                blocks.push(MovieBlock {
                    start_frame: start + 1,
                    end_frame: i,
                    start_position: cell_start_position(start, frame_count),
                    end_position: cell_end_position(i - 1, frame_count),
                    kind: kind.to_string(),
                    label: label.to_string(),
                });
            }
            start = i;
            current = next;
        }
    }

    blocks
}

fn build_ticks(frame_count: usize) -> Vec<MovieTick> {
    if frame_count <= 1 {
        return vec![MovieTick {
            frame: 1,
            label: "1".to_string(),
            position: frame_position(0, 1),
        }];
    }

    let raw_step = ((frame_count - 1) as f32 / 6.0).ceil().max(1.0) as usize;
    let step = nice_step(raw_step);
    let mut frames = vec![1usize];

    let mut frame = step;
    while frame <= frame_count {
        if frame != 1 {
            frames.push(frame);
        }
        frame += step;
    }

    if frames.last().copied() != Some(frame_count) {
        frames.push(frame_count);
    }

    frames
        .into_iter()
        .map(|frame| MovieTick {
            frame,
            label: frame.to_string(),
            position: frame_position(frame - 1, frame_count),
        })
        .collect()
}

fn nice_step(step: usize) -> usize {
    if step <= 1 {
        return 1;
    }

    let mut scale = 1usize;
    while scale * 10 < step {
        scale *= 10;
    }

    for unit in [1usize, 2, 3, 5, 10] {
        let candidate = unit * scale;
        if candidate >= step {
            return candidate;
        }
    }

    10 * scale
}

fn frame_position(frame_index: usize, frame_count: usize) -> f32 {
    if frame_count == 0 {
        return 0.5;
    }
    ((frame_index as f32 + 0.5) / frame_count as f32).clamp(0.0, 1.0)
}

fn cell_start_position(frame_index: usize, frame_count: usize) -> f32 {
    if frame_count == 0 {
        0.0
    } else {
        frame_index as f32 / frame_count as f32
    }
}

fn cell_end_position(frame_index: usize, frame_count: usize) -> f32 {
    if frame_count == 0 {
        1.0
    } else {
        ((frame_index + 1) as f32 / frame_count as f32).min(1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use patinae_scene::SceneView;

    #[test]
    fn empty_movie_has_default_lanes_and_one_state_block() {
        let movie = Movie::new();
        let timeline = MovieTimeline::from_movie(&movie);

        assert_eq!(timeline.frame_count, 1);
        assert_eq!(timeline.ticks.len(), 1);
        assert_eq!(timeline.lanes.len(), 4);
        assert_eq!(timeline.lane("mset").unwrap().blocks[0].label, "state 1");
    }

    #[test]
    fn trajectory_only_movie_uses_identity_state_blocks() {
        let mut movie = Movie::new();
        movie.set_n_object_states(3);

        let timeline = MovieTimeline::from_movie(&movie);
        let labels: Vec<_> = timeline
            .lane("mset")
            .unwrap()
            .blocks
            .iter()
            .map(|block| block.label.as_str())
            .collect();

        assert_eq!(timeline.frame_count, 3);
        assert_eq!(labels, vec!["state 1", "state 2", "state 3"]);
    }

    #[test]
    fn mset_runs_are_collapsed_into_blocks() {
        let mut movie = Movie::new();
        movie.set_from_spec(vec![1, 1, 2, 2, 2, 1]);

        let timeline = MovieTimeline::from_movie(&movie);
        let blocks = &timeline.lane("mset").unwrap().blocks;

        assert_eq!(blocks.len(), 3);
        assert_eq!((blocks[0].start_frame, blocks[0].end_frame), (1, 2));
        assert_eq!(blocks[0].label, "state 1");
        assert_eq!((blocks[1].start_frame, blocks[1].end_frame), (3, 5));
        assert_eq!(blocks[1].label, "state 2");
        assert_eq!((blocks[2].start_frame, blocks[2].end_frame), (6, 6));
    }

    #[test]
    fn camera_scene_and_object_keyframes_become_markers() {
        let mut movie = Movie::with_frames(4);
        movie.store_camera_keyframe(1, SceneView::default());
        movie.store_scene_keyframe(2, "binding", SceneView::default());
        movie.store_object_keyframe(3, "lig", Some(2), None);

        let timeline = MovieTimeline::from_movie(&movie);

        assert_eq!(timeline.lane("camera").unwrap().markers[0].frame, 2);
        assert!((timeline.lane("camera").unwrap().markers[0].position - 0.375).abs() < 0.001);
        assert_eq!(timeline.lane("scene").unwrap().markers[0].label, "binding");
        assert_eq!(timeline.lane("object").unwrap().markers[0].label, "lig");
        assert_eq!(
            timeline.lane("object").unwrap().blocks[0].label,
            "lig: state 2"
        );
    }

    #[test]
    fn ruler_ticks_use_compact_human_steps() {
        let mut movie = Movie::new();
        movie.set_frame_count(180);

        let timeline = MovieTimeline::from_movie(&movie);
        let labels: Vec<_> = timeline
            .ticks
            .iter()
            .map(|tick| tick.label.as_str())
            .collect();

        assert_eq!(labels, vec!["1", "30", "60", "90", "120", "150", "180"]);
    }
}
