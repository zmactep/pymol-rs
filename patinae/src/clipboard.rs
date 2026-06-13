//! System clipboard access for the native desktop host.

use copypasta::{ClipboardContext, ClipboardProvider};

/// Copy text into the system clipboard.
///
/// # Errors
///
/// Returns an error string when the platform clipboard cannot be opened or
/// written.
pub(crate) fn set_clipboard_text(text: &str) -> Result<(), String> {
    let mut clipboard =
        ClipboardContext::new().map_err(|err| format!("clipboard unavailable: {err}"))?;
    set_clipboard_text_with(&mut clipboard, text)
}

fn set_clipboard_text_with(
    clipboard: &mut impl ClipboardProvider,
    text: &str,
) -> Result<(), String> {
    clipboard
        .set_contents(text.to_owned())
        .map_err(|err| format!("clipboard write failed: {err}"))
}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use copypasta::ClipboardProvider;

    use super::*;

    type ClipboardResult<T> = Result<T, Box<dyn Error + Send + Sync + 'static>>;

    #[derive(Default)]
    struct FakeClipboard {
        contents: String,
        fail_write: bool,
    }

    impl ClipboardProvider for FakeClipboard {
        fn get_contents(&mut self) -> ClipboardResult<String> {
            Ok(self.contents.clone())
        }

        fn set_contents(&mut self, contents: String) -> ClipboardResult<()> {
            if self.fail_write {
                return Err(Box::<dyn Error + Send + Sync>::from("write denied"));
            }
            self.contents = contents;
            Ok(())
        }
    }

    #[test]
    fn set_clipboard_text_with_writes_text() {
        let mut clipboard = FakeClipboard::default();

        set_clipboard_text_with(&mut clipboard, ">obj:A\nACD").unwrap();

        assert_eq!(clipboard.contents, ">obj:A\nACD");
    }

    #[test]
    fn set_clipboard_text_with_reports_provider_error() {
        let mut clipboard = FakeClipboard {
            fail_write: true,
            ..FakeClipboard::default()
        };

        let err = set_clipboard_text_with(&mut clipboard, "ignored").unwrap_err();

        assert!(err.contains("clipboard write failed"));
        assert!(err.contains("write denied"));
    }
}
