//! Real filesystem answers for the questions the form asks.
//!
//! The design prototype faked all of this — it checked "does the string contain
//! a slash" and compared against a hard-coded list of three directories. These
//! are the real checks.

use std::path::{Path, PathBuf};

/// What the frontend needs to know about a path the user typed or picked.
#[derive(Debug, Default, Clone, serde::Serialize)]
pub struct PathInfo {
    pub exists: bool,
    pub is_dir: bool,
    pub is_file: bool,
    /// Lowercased extension without the dot ("raw", "mzml", "arrow", "pion").
    pub extension: String,
    /// Entries directly inside, when this is a directory.
    pub entry_count: usize,
    /// MS data files (.raw/.mzML/.arrow) directly inside, when a directory;
    /// 1 when the path is itself such a file.
    pub ms_file_count: usize,
    /// Counts per MS extension, so the UI can say "18 .raw files".
    pub raw_count: usize,
    pub mzml_count: usize,
    pub arrow_count: usize,
    /// True when this directory holds a config.json we could load.
    pub has_config_json: bool,
    /// True when this directory looks like an unpacked `.pion` spectral
    /// library. A `.pion` is a directory, not a file — the extension alone is
    /// not enough to tell a real library from an empty folder someone named
    /// `lib.pion`.
    pub is_pion_library: bool,
    /// Set when the path could not be read (permissions, broken symlink).
    pub error: Option<String>,
}

const MS_EXTENSIONS: [&str; 3] = ["raw", "mzml", "arrow"];

/// Files every `.pion` library contains. Checked by stem so the `.jld2`/`.jls`
/// serialization difference between library versions does not matter.
const PION_MARKERS: [&str; 2] = ["precursors_table", "detailed_fragments"];

/// Lowercased extension, looking through a trailing `.gz`.
///
/// Pioneer reads gzipped FASTA directly, so `proteins.fasta.gz` must report
/// "fasta" rather than "gz". `Path::extension` returns only the final component,
/// which made the validator reject every compressed FASTA even though the file
/// picker offered them and the search would have accepted them.
fn extension_of(p: &Path) -> String {
    let ext = |q: &Path| {
        q.extension()
            .and_then(|e| e.to_str())
            .map(|e| e.to_ascii_lowercase())
            .unwrap_or_default()
    };
    let last = ext(p);
    if last == "gz" {
        // `file_stem` of "proteins.fasta.gz" is "proteins.fasta"; take its
        // extension. Falls back to "gz" for a bare "archive.gz".
        let inner = p.file_stem().map(Path::new).map(ext).unwrap_or_default();
        if !inner.is_empty() {
            return inner;
        }
    }
    last
}

pub fn inspect(path: &str) -> PathInfo {
    let mut info = PathInfo::default();
    let trimmed = path.trim();
    if trimmed.is_empty() {
        return info;
    }
    let p = PathBuf::from(trimmed);

    let meta = match std::fs::metadata(&p) {
        Ok(m) => m,
        Err(e) => {
            // NotFound is the normal "still typing" case, not an error worth showing.
            if e.kind() != std::io::ErrorKind::NotFound {
                info.error = Some(e.to_string());
            }
            info.extension = extension_of(&p);
            return info;
        }
    };

    info.exists = true;
    info.is_dir = meta.is_dir();
    info.is_file = meta.is_file();
    info.extension = extension_of(&p);

    if info.is_file {
        if MS_EXTENSIONS.contains(&info.extension.as_str()) {
            info.ms_file_count = 1;
            match info.extension.as_str() {
                "raw" => info.raw_count = 1,
                "mzml" => info.mzml_count = 1,
                "arrow" => info.arrow_count = 1,
                _ => {}
            }
        }
        return info;
    }

    if info.is_dir {
        let mut markers_seen = 0usize;
        match std::fs::read_dir(&p) {
            Ok(entries) => {
                for entry in entries.flatten() {
                    info.entry_count += 1;
                    let name = entry.file_name();
                    let name = name.to_string_lossy();
                    if name.eq_ignore_ascii_case("config.json") {
                        info.has_config_json = true;
                    }
                    let stem = Path::new(name.as_ref())
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or_default();
                    if PION_MARKERS.contains(&stem) {
                        markers_seen += 1;
                    }
                    let ext = extension_of(Path::new(name.as_ref()));
                    match ext.as_str() {
                        "raw" => {
                            info.raw_count += 1;
                            info.ms_file_count += 1;
                        }
                        "mzml" => {
                            info.mzml_count += 1;
                            info.ms_file_count += 1;
                        }
                        "arrow" => {
                            info.arrow_count += 1;
                            info.ms_file_count += 1;
                        }
                        _ => {}
                    }
                }
            }
            Err(e) => info.error = Some(e.to_string()),
        }
        info.is_pion_library = markers_seen == PION_MARKERS.len();
    }

    info
}

/// Read a config.json, accepting either the file itself or a results folder
/// containing one — which is what the "Load previous run" dialog offers.
pub fn read_config(path: &str) -> Result<String, String> {
    let p = PathBuf::from(path.trim());
    let target = if p.is_dir() {
        let candidate = p.join("config.json");
        if !candidate.is_file() {
            return Err(format!("No config.json in {}", p.display()));
        }
        candidate
    } else {
        p
    };
    std::fs::read_to_string(&target).map_err(|e| format!("Could not read {}: {e}", target.display()))
}

#[cfg(test)]
mod extension_tests {
    use super::extension_of;
    use std::path::Path;

    #[test]
    fn looks_through_a_trailing_gz() {
        // The case Dennis hit: the picker offered it, the validator rejected it.
        assert_eq!(extension_of(Path::new("proteins.fasta.gz")), "fasta");
        assert_eq!(extension_of(Path::new("/a/b/UP000005640.fa.gz")), "fa");
        assert_eq!(extension_of(Path::new("x.FAA.GZ")), "faa");
    }

    #[test]
    fn leaves_uncompressed_paths_alone() {
        assert_eq!(extension_of(Path::new("proteins.fasta")), "fasta");
        assert_eq!(extension_of(Path::new("run.raw")), "raw");
        assert_eq!(extension_of(Path::new("lib.poin")), "poin");
        assert_eq!(extension_of(Path::new("noext")), "");
    }

    #[test]
    fn a_bare_gz_is_still_gz() {
        // Nothing to look through, so the trailing component stands.
        assert_eq!(extension_of(Path::new("archive.gz")), "gz");
    }

    #[test]
    fn a_dotfile_named_gz_has_no_extension() {
        // Rust treats a leading-dot name as a hidden file rather than an
        // extension, so ".gz" reports "" and never reaches the look-through.
        assert_eq!(extension_of(Path::new("/tmp/.gz")), "");
        assert_eq!(extension_of(Path::new(".fasta")), "");
    }
}

#[cfg(test)]
mod missing_path_tests {
    use super::inspect;

    /// A history entry can name a folder that has since been deleted, moved, or
    /// that lives on an unmounted volume. `inspect` must report that plainly
    /// rather than erroring: NotFound is the normal case here, not a fault.
    #[test]
    fn a_deleted_directory_reports_absent_without_an_error() {
        let dir = std::env::temp_dir().join("pioneer_gui_missing_test_lib.poin");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let before = inspect(dir.to_str().unwrap());
        assert!(before.exists, "fixture should exist before removal");
        assert!(before.is_dir);

        std::fs::remove_dir_all(&dir).unwrap();

        let after = inspect(dir.to_str().unwrap());
        assert!(!after.exists, "a removed directory must report exists = false");
        assert!(!after.is_dir);
        assert!(
            after.error.is_none(),
            "NotFound is expected, not an error to surface: {:?}",
            after.error
        );
        // The extension still parses, so the UI can keep describing what it was.
        assert_eq!(after.extension, "poin");
    }

    #[test]
    fn a_library_folder_stops_looking_like_one_when_emptied() {
        let dir = std::env::temp_dir().join("pioneer_gui_emptied_lib.poin");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("precursors_table.arrow"), b"x").unwrap();
        std::fs::write(dir.join("detailed_fragments.jls"), b"x").unwrap();
        assert!(inspect(dir.to_str().unwrap()).is_pion_library);

        // A folder that still exists but has lost its contents is the nastier
        // case: the path resolves, so only the marker check catches it.
        std::fs::remove_file(dir.join("precursors_table.arrow")).unwrap();
        let after = inspect(dir.to_str().unwrap());
        assert!(after.exists);
        assert!(!after.is_pion_library, "an emptied folder must not pass as a library");

        std::fs::remove_dir_all(&dir).unwrap();
    }
}
