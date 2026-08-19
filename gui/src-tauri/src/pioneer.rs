//! Locating the packaged Pioneer distribution and running its executables.
//!
//! The distribution layout we expect (as produced by PackageCompiler, plus the
//! .NET converter added at packaging time):
//!
//! ```text
//! <pioneer_home>/
//!     pioneer                 wrapper script: maps a subcommand to bin/<exe> and execs it
//!     bin/SearchDIA
//!     bin/BuildSpecLib
//!     bin/GetSearchParams
//!     bin/GetBuildLibParams
//!     bin/convertMzML
//!     bin/PioneerConverter    separate .NET 8 repo
//!     lib/
//! ```

use std::path::{Path, PathBuf};

/// The Pioneer subcommands the GUI can drive, and the executable each maps to.
///
/// We invoke `bin/<exe>` directly rather than going through the `pioneer`
/// wrapper: the wrapper hard-codes `JULIA_NUM_THREADS=auto`, and the GUI needs
/// to pass the user's chosen thread count instead. Bypassing it means *we* are
/// responsible for setting `JULIA_NUM_THREADS` — a bare `bin/SearchDIA` with no
/// such variable silently runs single-threaded (~3.7x slower).
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Command {
    SearchDia,
    BuildSpecLib,
    DownloadSpecLib,
    ConvertRaw,
}

impl Command {
    pub fn exe_name(self) -> &'static str {
        match self {
            Command::SearchDia => "SearchDIA",
            Command::BuildSpecLib => "BuildSpecLib",
            Command::DownloadSpecLib => "DownloadSpecLib",
            Command::ConvertRaw => "PioneerConverter",
        }
    }

    /// The subcommand the `pioneer` wrapper uses, for the fallback path.
    pub fn subcommand(self) -> &'static str {
        match self {
            Command::SearchDia => "search",
            Command::BuildSpecLib => "predict",
            Command::DownloadSpecLib => "download",
            Command::ConvertRaw => "convert-raw",
        }
    }

    /// Julia-based executables need JULIA_NUM_THREADS; the .NET converter does not.
    pub fn is_julia(self) -> bool {
        !matches!(self, Command::ConvertRaw)
    }
}

/// Where a distribution was found, and what it contains.
#[derive(Debug, Clone, serde::Serialize)]
pub struct PioneerInfo {
    /// Root of the distribution.
    pub home: String,
    /// How we found it — shown in the UI so a wrong pickup is diagnosable.
    pub source: String,
    /// Executables present under `bin/`, by name.
    pub executables: Vec<String>,
    /// True when the `pioneer` wrapper script is present.
    pub has_wrapper: bool,
    /// Contents of the distribution's `VERSION` file, when it has one. The CLI
    /// tools and the converter ship together, so this one string versions
    /// both. None for a distribution built before the file existed.
    pub version: Option<String>,
}

/// Where the platform installer puts the CLI that belongs to this GUI build.
///
/// Release builds set `PIONEER_INSTALL_KEY` while compiling the GUI. Keeping
/// the key in the binary means `Pioneer 2.1.0.app` always launches Pioneer
/// 2.1.0 even when several versions are installed and another version owns the
/// unqualified `pioneer` command.
fn installed_home() -> Option<PathBuf> {
    let install_key = option_env!("PIONEER_INSTALL_KEY")
        .map(str::trim)
        .filter(|value| !value.is_empty());

    #[cfg(windows)]
    {
        let pf = std::env::var_os("ProgramFiles")
            .unwrap_or_else(|| std::ffi::OsString::from("C:\\Program Files"));
        let root = PathBuf::from(pf).join("Pioneer");
        Some(install_key.map_or(root.clone(), |key| root.join(key)))
    }
    #[cfg(target_os = "macos")]
    {
        let root = PathBuf::from("/usr/local/lib/pioneer");
        // Installer builds embed their version key and must remain pinned to
        // that distribution. An unversioned `tauri dev` build instead follows
        // the CLI's active-version link, just like /usr/local/bin/pioneer.
        Some(install_key.map_or_else(|| root.join("current"), |key| root.join(key)))
    }
    #[cfg(all(unix, not(target_os = "macos")))]
    {
        let root = PathBuf::from("/opt/pioneer");
        Some(install_key.map_or(root.clone(), |key| root.join(key)))
    }
}

/// Resolve the distribution root, in priority order:
///
/// 1. `PIONEER_HOME` — the dev-loop override, and the escape hatch when the
///    installed copy is missing or wrong.
/// 2. The platform install location the Pioneer installer writes to.
/// 3. `<resources>/pioneer` — only present if someone builds a deliberately
///    self-contained bundle; the shipped installer does not populate it.
///
/// Returns the first candidate that looks like a real distribution.
pub fn resolve_home(resource_dir: Option<&Path>) -> Result<PioneerInfo, String> {
    resolve_home_with_installed(resource_dir, installed_home())
}

fn resolve_home_with_installed(
    resource_dir: Option<&Path>,
    installed: Option<PathBuf>,
) -> Result<PioneerInfo, String> {
    let mut tried: Vec<String> = Vec::new();

    if let Some(env_home) = std::env::var_os("PIONEER_HOME") {
        let p = PathBuf::from(env_home);
        match inspect(&p, "PIONEER_HOME") {
            Some(info) => return Ok(info),
            None => tried.push(format!("PIONEER_HOME={}", p.display())),
        }
    }

    if let Some(p) = installed {
        match inspect(&p, "installed with Pioneer") {
            Some(info) => return Ok(info),
            None => tried.push(format!("installed: {}", p.display())),
        }
    }

    if let Some(res) = resource_dir {
        let p = res.join("pioneer");
        match inspect(&p, "bundled with the app") {
            Some(info) => return Ok(info),
            None => tried.push(format!("bundled: {}", p.display())),
        }
    }

    Err(format!(
        "No Pioneer distribution found. Looked in: {}. \
         Install Pioneer, or set PIONEER_HOME to an unpacked distribution.",
        if tried.is_empty() { "nowhere — no candidates configured".into() } else { tried.join("; ") }
    ))
}

/// Return distribution info if `home` contains at least one known executable.
fn inspect(home: &Path, source: &str) -> Option<PioneerInfo> {
    if !home.is_dir() {
        return None;
    }
    let bin = home.join("bin");
    let mut executables = Vec::new();
    for cmd in [Command::SearchDia, Command::BuildSpecLib, Command::DownloadSpecLib,
                Command::ConvertRaw] {
        if exe_path(&bin, cmd.exe_name()).is_some() {
            executables.push(cmd.exe_name().to_string());
        }
    }
    let has_wrapper = exe_path(home, "pioneer").is_some();

    if executables.is_empty() && !has_wrapper {
        return None;
    }
    let version = std::fs::read_to_string(home.join("VERSION"))
        .ok()
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty());

    Some(PioneerInfo {
        home: home.display().to_string(),
        source: source.to_string(),
        executables,
        has_wrapper,
        version,
    })
}

/// Look for `name` (and `name.exe` on Windows) inside `dir`.
fn exe_path(dir: &Path, name: &str) -> Option<PathBuf> {
    let direct = dir.join(name);
    if direct.is_file() {
        return Some(direct);
    }
    let with_ext = dir.join(format!("{name}.exe"));
    if with_ext.is_file() {
        return Some(with_ext);
    }
    None
}

/// Resolve the executable to spawn for `cmd`.
///
/// Prefers `bin/<exe>` so we control `JULIA_NUM_THREADS`; falls back to the
/// `pioneer` wrapper (which sets `JULIA_NUM_THREADS=auto` itself, ignoring the
/// UI's thread count — the caller is told so it can warn).
pub struct Resolved {
    pub program: PathBuf,
    /// Arguments that must precede the params file (the wrapper's subcommand).
    pub leading_args: Vec<String>,
    /// True when we fell back to the wrapper and cannot control thread count.
    pub via_wrapper: bool,
}

pub fn resolve_command(home: &Path, cmd: Command) -> Result<Resolved, String> {
    if let Some(p) = exe_path(&home.join("bin"), cmd.exe_name()) {
        return Ok(Resolved { program: p, leading_args: Vec::new(), via_wrapper: false });
    }
    if let Some(p) = exe_path(home, "pioneer") {
        return Ok(Resolved {
            program: p,
            leading_args: vec![cmd.subcommand().to_string()],
            via_wrapper: true,
        });
    }
    Err(format!(
        "{} not found in {} (looked for bin/{} and the pioneer wrapper)",
        cmd.exe_name(),
        home.display(),
        cmd.exe_name()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A development test build has no embedded install key and therefore uses
    /// the platform's unversioned/active location. Installer builds append
    /// their compile-time key.
    #[test]
    fn installed_home_matches_the_installer_layout() {
        let home = installed_home().expect("a platform default");
        #[cfg(target_os = "macos")]
        {
            assert!(home.starts_with("/usr/local/lib/pioneer"));
            if option_env!("PIONEER_INSTALL_KEY").is_none() {
                assert!(
                    home.ends_with("current"),
                    "dev build should follow {home:?}"
                );
            }
        }
        #[cfg(all(unix, not(target_os = "macos")))]
        assert!(home.starts_with("/opt/pioneer"));
        #[cfg(windows)]
        assert!(home.to_string_lossy().contains("Pioneer"));
        if let Some(key) = option_env!("PIONEER_INSTALL_KEY") {
            assert!(home.ends_with(key), "installer key missing from {home:?}");
        }
    }

    /// With nothing configured, the error must name every place we looked —
    /// a silent "not found" is what makes a bad install hard to diagnose.
    #[test]
    fn error_lists_every_candidate() {
        // SAFETY: single-threaded test, restored immediately.
        unsafe { std::env::remove_var("PIONEER_HOME") };
        let missing = std::env::temp_dir().join("pioneer-gui-missing-install-for-test");
        let err = resolve_home_with_installed(None, Some(missing)).unwrap_err();
        assert!(err.contains("installed:"), "should report the install path: {err}");
        assert!(err.contains("PIONEER_HOME"), "should mention the override: {err}");
    }
}
