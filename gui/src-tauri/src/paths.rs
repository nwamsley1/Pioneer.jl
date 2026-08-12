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

/// What a `.poin` library says about itself.
///
/// Read from the library's own `config.json` — the BuildSpecLib config that
/// produced it — rather than parsed out of its filename. A name breaks the
/// moment a library is renamed or built elsewhere, and cannot report anything
/// the name does not already say.
///
/// Every field is optional in practice: libraries built before a key existed
/// simply lack it, so each is reported as empty rather than guessed at. The
/// oldest libraries have no `prediction_model` at all.
#[derive(Debug, Default, Clone, serde::Serialize)]
pub struct LibraryInfo {
    /// False when the path is not a library, or has no readable config.
    pub is_library: bool,
    pub prediction_model: String,
    /// "7–30" — empty when the digest parameters are absent.
    pub length_range: String,
    pub charge_range: String,
    pub missed_cleavages: String,
    pub max_var_mods: String,
    /// "Unimod:35 on M" per entry.
    pub fixed_mods: Vec<String>,
    pub variable_mods: Vec<String>,
    pub fastas: Vec<String>,
    pub include_contaminants: bool,
    pub has_decoys: bool,
    pub nce: String,
    /// Set when the library exists but its config could not be read.
    pub error: Option<String>,
}

/// Format one `{pattern, mass, name}` block into "name on pattern" strings.
///
/// The three arrays are parallel; a short or missing one is tolerated rather
/// than panicking, since this is describing someone else's file.
fn mod_list(v: Option<&serde_json::Value>) -> Vec<String> {
    let Some(obj) = v else { return Vec::new() };
    let names = obj.get("name").and_then(|x| x.as_array());
    let patterns = obj.get("pattern").and_then(|x| x.as_array());
    let Some(names) = names else { return Vec::new() };
    names
        .iter()
        .enumerate()
        .map(|(i, n)| {
            let name = n.as_str().unwrap_or("?");
            match patterns.and_then(|p| p.get(i)).and_then(|p| p.as_str()) {
                Some(pat) => format!("{name} on {pat}"),
                None => name.to_string(),
            }
        })
        .collect()
}

fn num(v: Option<&serde_json::Value>) -> String {
    match v {
        Some(x) if x.is_number() => x.to_string(),
        _ => String::new(),
    }
}

pub fn library_info(path: &str) -> LibraryInfo {
    let mut out = LibraryInfo::default();
    let p = PathBuf::from(path.trim());
    if !inspect(path).is_pion_library {
        return out;
    }
    out.is_library = true;

    let text = match std::fs::read_to_string(p.join("config.json")) {
        Ok(t) => t,
        Err(e) => {
            out.error = Some(format!("No readable config.json: {e}"));
            return out;
        }
    };
    let cfg: serde_json::Value = match serde_json::from_str(&text) {
        Ok(v) => v,
        Err(e) => {
            out.error = Some(format!("config.json is not valid JSON: {e}"));
            return out;
        }
    };

    if let Some(m) = cfg
        .get("library_params")
        .and_then(|l| l.get("prediction_model"))
        .and_then(|m| m.as_str())
    {
        out.prediction_model = m.to_string();
    }

    if let Some(d) = cfg.get("fasta_digest_params") {
        let (lo, hi) = (num(d.get("min_length")), num(d.get("max_length")));
        if !lo.is_empty() && !hi.is_empty() {
            out.length_range = format!("{lo}\u{2013}{hi}");
        }
        let (clo, chi) = (num(d.get("min_charge")), num(d.get("max_charge")));
        if !clo.is_empty() && !chi.is_empty() {
            out.charge_range = format!("{clo}\u{2013}{chi}");
        }
        out.missed_cleavages = num(d.get("missed_cleavages"));
        out.max_var_mods = num(d.get("max_var_mods"));
        out.has_decoys = d.get("add_decoys").and_then(|x| x.as_bool()).unwrap_or(false);
    }

    out.nce = num(cfg.get("nce_params").and_then(|n| n.get("nce")));
    out.fixed_mods = mod_list(cfg.get("fixed_mods"));
    out.variable_mods = mod_list(cfg.get("variable_mods"));
    out.include_contaminants = cfg
        .get("include_contaminants")
        .and_then(|x| x.as_bool())
        .unwrap_or(false);
    if let Some(names) = cfg.get("fasta_names").and_then(|x| x.as_array()) {
        out.fastas = names
            .iter()
            .filter_map(|n| n.as_str().map(str::to_string))
            .collect();
    }
    out
}

#[cfg(test)]
mod library_info_tests {
    use super::library_info;

    #[test]
    fn a_non_library_reports_nothing() {
        let d = std::env::temp_dir().join("pioneer_libinfo_not_a_lib");
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        let info = library_info(d.to_str().unwrap());
        assert!(!info.is_library);
        assert!(info.error.is_none(), "not-a-library is not an error");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn a_library_without_a_config_says_so() {
        // The markers make it a library; the missing config is the fault.
        let d = std::env::temp_dir().join("pioneer_libinfo_no_config.poin");
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        std::fs::write(d.join("precursors_table.arrow"), b"x").unwrap();
        std::fs::write(d.join("detailed_fragments.jls"), b"x").unwrap();
        let info = library_info(d.to_str().unwrap());
        assert!(info.is_library);
        assert!(info.error.is_some(), "a library with no config should report why");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn reads_what_the_config_records_and_leaves_the_rest_empty() {
        let d = std::env::temp_dir().join("pioneer_libinfo_real.poin");
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        std::fs::write(d.join("precursors_table.arrow"), b"x").unwrap();
        std::fs::write(d.join("detailed_fragments.jls"), b"x").unwrap();
        // No prediction_model: the shape older libraries actually have.
        std::fs::write(
            d.join("config.json"),
            br#"{
              "library_params": {"frag_mz_min": 150.0},
              "fasta_digest_params": {"min_length": 7, "max_length": 30,
                                      "min_charge": 2, "max_charge": 4,
                                      "missed_cleavages": 1, "max_var_mods": 1,
                                      "add_decoys": true},
              "nce_params": {"nce": 26.0},
              "variable_mods": {"pattern": ["M"], "name": ["Unimod:35"], "mass": [15.99491]},
              "fixed_mods": {"pattern": ["C"], "name": ["Unimod:4"], "mass": [57.021464]},
              "fasta_names": ["ECOLI", "CONTAM"],
              "include_contaminants": true
            }"#,
        )
        .unwrap();

        let i = library_info(d.to_str().unwrap());
        assert!(i.is_library);
        assert!(i.error.is_none());
        assert_eq!(i.length_range, "7\u{2013}30");
        assert_eq!(i.charge_range, "2\u{2013}4");
        assert_eq!(i.missed_cleavages, "1");
        assert_eq!(i.variable_mods, vec!["Unimod:35 on M"]);
        assert_eq!(i.fixed_mods, vec!["Unimod:4 on C"]);
        assert_eq!(i.fastas, vec!["ECOLI", "CONTAM"]);
        assert!(i.include_contaminants);
        assert!(i.has_decoys);
        assert_eq!(i.nce, "26.0");
        // Absent rather than guessed: the oldest libraries do not record it.
        assert_eq!(i.prediction_model, "");

        let _ = std::fs::remove_dir_all(&d);
    }
}

#[cfg(test)]
mod real_library_tests {
    use super::library_info;

    /// Reads the actual libraries on this machine when they are present, so the
    /// parser is checked against real files rather than only a fixture.
    #[test]
    fn reads_the_libraries_on_disk() {
        let home = std::env::var("HOME").unwrap_or_default();
        for (path, expect_model) in [
            (format!("{home}/Desktop/PioneerTest/ecoli.poin"), ""),
            (
                format!("{home}/Projects/Pioneer.jl/data/precompile/ecoli_small_prosit.poin"),
                "prosit_2020_hcd",
            ),
        ] {
            if !std::path::Path::new(&path).is_dir() {
                continue;
            }
            let i = library_info(&path);
            assert!(i.is_library, "{path} should read as a library");
            assert!(i.error.is_none(), "{path}: {:?}", i.error);
            assert_eq!(i.prediction_model, expect_model, "model for {path}");
            assert!(!i.length_range.is_empty(), "{path} should record a length range");
            assert!(!i.fixed_mods.is_empty(), "{path} should record fixed mods");
            println!(
                "  {path}\n    model={:?} digest={} charge={} mods={:?} fasta={:?}",
                i.prediction_model, i.length_range, i.charge_range, i.variable_mods, i.fastas
            );
        }
    }
}
