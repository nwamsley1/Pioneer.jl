//! macOS package uninstall support.
//!
//! The package installs a private, version-specific helper under the CLI root.
//! Development builds and non-macOS packages do not expose an uninstall action.

use std::path::PathBuf;

#[derive(Debug, Clone, serde::Serialize)]
pub struct Info {
    pub available: bool,
    pub version: String,
    pub install_root: String,
    pub app_path: String,
}

fn install_key() -> Option<&'static str> {
    option_env!("PIONEER_INSTALL_KEY")
        .map(str::trim)
        .filter(|value| valid_install_key(value))
}

fn valid_install_key(value: &str) -> bool {
    !value.is_empty()
        && value
            .chars()
            .next()
            .is_some_and(|c| c.is_ascii_alphanumeric())
        && value
            .chars()
            .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '.' || c == '-')
}

#[cfg(target_os = "macos")]
fn target() -> Option<(Info, PathBuf)> {
    let key = install_key()?;
    let root = PathBuf::from("/usr/local/lib/pioneer").join(key);
    let helper = root.join("libexec/uninstall");
    let version = std::fs::read_to_string(root.join("VERSION"))
        .ok()?
        .trim()
        .to_string();
    if version.is_empty()
        || !version
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '+' | '-'))
    {
        return None;
    }
    let app_path = PathBuf::from("/Applications").join(format!("Pioneer {version}.app"));
    let available = root.is_dir() && helper.is_file() && app_path.is_dir();
    Some((
        Info {
            available,
            version,
            install_root: root.display().to_string(),
            app_path: app_path.display().to_string(),
        },
        helper,
    ))
}

pub fn info() -> Info {
    #[cfg(target_os = "macos")]
    if let Some((info, _)) = target() {
        return info;
    }

    Info {
        available: false,
        version: String::new(),
        install_root: String::new(),
        app_path: String::new(),
    }
}

#[cfg(target_os = "macos")]
pub async fn run() -> Result<(), String> {
    let (info, helper) = target().ok_or_else(|| {
        "This copy of Pioneer was not installed by a macOS Pioneer package.".to_string()
    })?;
    if !info.available {
        return Err("The matching Pioneer installation is incomplete or has moved.".to_string());
    }

    tauri::async_runtime::spawn_blocking(move || {
        let output = std::process::Command::new("/usr/bin/osascript")
            .args([
                "-e",
                "on run argv",
                "-e",
                "do shell script (quoted form of (item 1 of argv)) with administrator privileges",
                "-e",
                "end run",
                "--",
            ])
            .arg(&helper)
            .output()
            .map_err(|e| format!("Could not open the macOS administrator prompt: {e}"))?;

        if output.status.success() {
            return Ok(());
        }
        let detail = String::from_utf8_lossy(&output.stderr).trim().to_string();
        if detail.contains("User canceled") || detail.contains("(-128)") {
            Err("Uninstall was cancelled.".to_string())
        } else if detail.is_empty() {
            Err(format!("Uninstall failed with status {}.", output.status))
        } else {
            Err(format!("Uninstall failed: {detail}"))
        }
    })
    .await
    .map_err(|e| format!("Uninstall task failed: {e}"))?
}

#[cfg(not(target_os = "macos"))]
pub async fn run() -> Result<(), String> {
    Err("Uninstall from the GUI is only available in the packaged macOS app.".to_string())
}

#[cfg(test)]
mod tests {
    use super::valid_install_key;

    #[test]
    fn accepts_only_normalized_install_keys() {
        assert!(valid_install_key("2.1.0"));
        assert!(valid_install_key("2.1.0-rc1"));
        assert!(!valid_install_key(""));
        assert!(!valid_install_key("../2.1.0"));
        assert!(!valid_install_key("V2.1.0"));
        assert!(!valid_install_key("2.1.0 rc1"));
    }
}
