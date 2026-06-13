//! Python syntax highlighting for the scripting panel.

use patinae_plugin::prelude::{PanelTextHighlight, PanelTextStyle};
use tree_sitter_highlight::{HighlightConfiguration, HighlightEvent, Highlighter};

const HIGHLIGHT_NAMES: &[&str] = &[
    "attribute",
    "comment",
    "constant",
    "constant.builtin",
    "constructor",
    "function",
    "function.builtin",
    "keyword",
    "number",
    "operator",
    "property",
    "punctuation",
    "punctuation.bracket",
    "punctuation.delimiter",
    "string",
    "string.special",
    "type",
    "type.builtin",
    "variable.builtin",
];

#[derive(Debug, Default)]
pub(crate) struct PythonHighlightCache {
    source: String,
    spans: Vec<PanelTextHighlight>,
}

impl PythonHighlightCache {
    pub(crate) fn highlights_for(&mut self, source: &str) -> Vec<PanelTextHighlight> {
        if self.source != source {
            self.source = source.to_string();
            self.spans = highlight_python(source);
        }
        self.spans.clone()
    }
}

pub(crate) fn highlight_python(source: &str) -> Vec<PanelTextHighlight> {
    if source.is_empty() {
        return Vec::new();
    }

    let Ok(mut config) = HighlightConfiguration::new(
        tree_sitter_python::LANGUAGE.into(),
        "python",
        tree_sitter_python::HIGHLIGHTS_QUERY,
        "",
        "",
    ) else {
        return Vec::new();
    };
    config.configure(HIGHLIGHT_NAMES);

    let mut highlighter = Highlighter::new();
    let Ok(events) = highlighter.highlight(&config, source.as_bytes(), None, |_| None) else {
        return Vec::new();
    };

    let mut active = Vec::new();
    let mut spans = Vec::new();
    for event in events {
        match event {
            Ok(HighlightEvent::HighlightStart(highlight)) => {
                active.push(style_for_highlight(highlight.0));
            }
            Ok(HighlightEvent::HighlightEnd) => {
                active.pop();
            }
            Ok(HighlightEvent::Source { start, end }) => {
                let Some(style) = active.iter().rev().find_map(|style| *style) else {
                    continue;
                };
                if start < end && end <= source.len() {
                    spans.push(PanelTextHighlight::new(start, end, style));
                }
            }
            Err(_) => return Vec::new(),
        }
    }

    spans
}

fn style_for_highlight(index: usize) -> Option<PanelTextStyle> {
    match HIGHLIGHT_NAMES.get(index).copied()? {
        "comment" => Some(PanelTextStyle::Comment),
        "constant" | "constant.builtin" => Some(PanelTextStyle::Constant),
        "constructor" | "type" | "type.builtin" => Some(PanelTextStyle::Type),
        "function" => Some(PanelTextStyle::Function),
        "function.builtin" | "variable.builtin" => Some(PanelTextStyle::Builtin),
        "keyword" => Some(PanelTextStyle::Keyword),
        "number" => Some(PanelTextStyle::Number),
        "operator" => Some(PanelTextStyle::Operator),
        "punctuation" | "punctuation.bracket" | "punctuation.delimiter" => {
            Some(PanelTextStyle::Punctuation)
        }
        "string" | "string.special" => Some(PanelTextStyle::String),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn has_style(source: &str, spans: &[PanelTextHighlight], style: PanelTextStyle) -> bool {
        spans.iter().any(|span| {
            span.style == style
                && span.start < span.end
                && span.end <= source.len()
                && source.is_char_boundary(span.start)
                && source.is_char_boundary(span.end)
        })
    }

    #[test]
    fn highlights_python_constructs() {
        let source = "def f(x):\n    return \"hi\" # note\n";
        let spans = highlight_python(source);

        assert!(has_style(source, &spans, PanelTextStyle::Keyword));
        assert!(has_style(source, &spans, PanelTextStyle::Function));
        assert!(has_style(source, &spans, PanelTextStyle::String));
        assert!(has_style(source, &spans, PanelTextStyle::Comment));
    }

    #[test]
    fn spans_stay_inside_utf8_boundaries() {
        let source = "name = 'привет'\n# note\n";
        for span in highlight_python(source) {
            assert!(span.start < span.end);
            assert!(span.end <= source.len());
            assert!(source.is_char_boundary(span.start));
            assert!(source.is_char_boundary(span.end));
        }
    }

    #[test]
    fn incomplete_python_returns_safe_result() {
        let source = "def broken(:\n    print('x'";
        for span in highlight_python(source) {
            assert!(span.start < span.end);
            assert!(span.end <= source.len());
        }
    }
}
