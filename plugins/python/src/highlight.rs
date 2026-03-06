//! Python Syntax Highlighting
//!
//! Uses syntect to tokenize Python code and produce an egui `LayoutJob`
//! with colored text spans for use with `TextEdit::layouter`.

use egui::text::LayoutJob;
use egui::{Color32, FontId, TextFormat};
use syntect::easy::HighlightLines;
use syntect::highlighting::{self, ThemeSet};
use syntect::parsing::SyntaxSet;

/// Cached syntax highlighting state for Python.
pub struct PythonHighlighter {
    syntax_set: SyntaxSet,
    theme: highlighting::Theme,
}

impl PythonHighlighter {
    pub fn new() -> Self {
        let syntax_set = SyntaxSet::load_defaults_newlines();
        let theme_set = ThemeSet::load_defaults();
        let theme = theme_set.themes["base16-mocha.dark"].clone();
        Self { syntax_set, theme }
    }

    /// Produce a colored `LayoutJob` from Python source code.
    pub fn highlight(&self, text: &str, font_id: FontId, wrap_width: f32) -> LayoutJob {
        let mut job = LayoutJob::default();
        job.wrap.max_width = wrap_width;

        let syntax = self
            .syntax_set
            .find_syntax_by_extension("py")
            .unwrap_or_else(|| self.syntax_set.find_syntax_plain_text());

        let mut highlighter = HighlightLines::new(syntax, &self.theme);

        for line in text.split_inclusive('\n') {
            let ranges = highlighter
                .highlight_line(line, &self.syntax_set)
                .unwrap_or_default();

            for (style, piece) in ranges {
                let fg = style.foreground;
                let color = Color32::from_rgb(fg.r, fg.g, fg.b);
                job.append(
                    piece,
                    0.0,
                    TextFormat {
                        font_id: font_id.clone(),
                        color,
                        ..Default::default()
                    },
                );
            }
        }

        // Handle empty text — append an empty span so the layout has valid metrics
        if text.is_empty() {
            job.append(
                "",
                0.0,
                TextFormat {
                    font_id,
                    ..Default::default()
                },
            );
        }

        job
    }
}
