use std::collections::HashMap;
use std::io;
use std::path::{Path, PathBuf};
use std::rc::Rc;

use slint::{ComponentHandle, Model, ModelRc, VecModel};

use patinae_settings::recent_files::RecentFiles;

use crate::recent_thumbnails::{extension_badge, RecentThumbnailCache, ThumbnailWritePaths};
use crate::{AppWindow, MenuState, RecentFileItem};

#[derive(Debug, Clone, PartialEq, Eq)]
struct RecentFileDisplay {
    title: String,
    subtitle: String,
    path: String,
    extension: String,
}

pub(crate) struct RecentFilesBridge {
    storage_path: PathBuf,
    store: RecentFiles,
    thumbnails: RecentThumbnailCache,
    model: Rc<VecModel<RecentFileItem>>,
}

impl RecentFilesBridge {
    pub(crate) fn new(storage_path: PathBuf, thumbnail_dir: PathBuf) -> Self {
        let store = RecentFiles::load(&storage_path);
        let thumbnails = RecentThumbnailCache::new(thumbnail_dir);
        let model = Rc::new(VecModel::from(slint_items(store.paths(), &thumbnails)));
        Self {
            storage_path,
            store,
            thumbnails,
            model,
        }
    }

    pub(crate) fn attach(&self, window: &AppWindow) {
        let menu = window.global::<MenuState>();
        menu.set_recent_files(ModelRc::from(self.model.clone()));
        menu.set_has_recent_files(self.model.row_count() > 0);
    }

    pub(crate) fn record(&mut self, path: PathBuf, window: &AppWindow) -> io::Result<()> {
        let previous = self.store.clone();
        self.store.record(path);
        if let Err(err) = self.store.save(&self.storage_path) {
            self.store = previous;
            return Err(err);
        }
        self.refresh(window);
        Ok(())
    }

    pub(crate) fn remove(&mut self, path: &Path, window: &AppWindow) -> io::Result<bool> {
        let previous = self.store.clone();
        let removed = self.store.remove(path);
        if removed {
            if let Err(err) = self.store.save(&self.storage_path) {
                self.store = previous;
                return Err(err);
            }
            self.refresh(window);
        }
        Ok(removed)
    }

    pub(crate) fn refresh(&self, window: &AppWindow) {
        replace_model(
            &self.model,
            slint_items(self.store.paths(), &self.thumbnails),
        );
        window
            .global::<MenuState>()
            .set_has_recent_files(self.model.row_count() > 0);
    }

    pub(crate) fn thumbnail_write_paths(&self, path: &Path) -> io::Result<ThumbnailWritePaths> {
        self.thumbnails.write_paths(path)
    }
}

fn slint_items(paths: &[PathBuf], thumbnails: &RecentThumbnailCache) -> Vec<RecentFileItem> {
    display_items(paths)
        .into_iter()
        .map(|item| {
            let thumbnail = thumbnails.load(Path::new(&item.path));
            let has_thumbnail = thumbnail.is_some();
            RecentFileItem {
                thumbnail: thumbnail.unwrap_or_default(),
                has_thumbnail,
                title: item.title.into(),
                subtitle: item.subtitle.into(),
                path: item.path.into(),
                extension: item.extension.into(),
            }
        })
        .collect()
}

fn display_items(paths: &[PathBuf]) -> Vec<RecentFileDisplay> {
    let mut basename_counts = HashMap::<String, usize>::new();
    for path in paths {
        let basename = display_basename(path);
        *basename_counts.entry(basename).or_default() += 1;
    }

    paths
        .iter()
        .map(|path| {
            let basename = display_basename(path);
            let parent = display_parent(path);
            let title = if basename_counts.get(&basename).copied().unwrap_or_default() > 1 {
                format!("{basename} - {parent}")
            } else {
                basename
            };
            RecentFileDisplay {
                title,
                subtitle: parent,
                path: path.to_string_lossy().into_owned(),
                extension: extension_badge(path),
            }
        })
        .collect()
}

fn display_basename(path: &Path) -> String {
    path.file_name()
        .and_then(|name| name.to_str())
        .map(str::to_string)
        .unwrap_or_else(|| path.display().to_string())
}

fn display_parent(path: &Path) -> String {
    path.parent()
        .and_then(|parent| {
            parent
                .file_name()
                .and_then(|name| name.to_str())
                .map(str::to_string)
                .or_else(|| Some(parent.display().to_string()).filter(|s| !s.is_empty()))
        })
        .unwrap_or_default()
}

fn replace_model<T: Clone + 'static>(model: &Rc<VecModel<T>>, rows: Vec<T>) {
    model.set_vec(rows);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_items_use_basename_for_unique_files() {
        let items = display_items(&[
            PathBuf::from("/tmp/alpha/1fsd.pdb"),
            PathBuf::from("/tmp/beta/1ubq.cif"),
        ]);

        assert_eq!(items[0].title, "1fsd.pdb");
        assert_eq!(items[0].subtitle, "alpha");
        assert_eq!(items[1].title, "1ubq.cif");
        assert_eq!(items[1].subtitle, "beta");
    }

    #[test]
    fn display_items_add_parent_context_for_duplicate_basenames() {
        let items = display_items(&[
            PathBuf::from("/tmp/alpha/model.pdb"),
            PathBuf::from("/tmp/beta/model.pdb"),
        ]);

        assert_eq!(items[0].title, "model.pdb - alpha");
        assert_eq!(items[1].title, "model.pdb - beta");
    }

    #[test]
    fn display_items_include_extension_badges_for_non_default_formats() {
        let items = display_items(&[
            PathBuf::from("/tmp/alpha/model.pdb"),
            PathBuf::from("/tmp/beta/session.pse"),
            PathBuf::from("/tmp/gamma/ligand.mol2"),
        ]);

        assert_eq!(items[0].extension, "");
        assert_eq!(items[1].extension, ".pse");
        assert_eq!(items[2].extension, ".mol2");
    }

    #[test]
    fn slint_items_use_fallback_when_thumbnail_is_missing() {
        let cache = RecentThumbnailCache::new(PathBuf::from("/tmp/patinae-missing-thumbnails"));
        let rows = slint_items(&[PathBuf::from("/tmp/alpha/session.pse")], &cache);

        assert_eq!(rows[0].title.as_str(), "session.pse");
        assert_eq!(rows[0].extension.as_str(), ".pse");
        assert!(!rows[0].has_thumbnail);
    }
}
