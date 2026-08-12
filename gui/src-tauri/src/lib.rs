mod paths;
mod pioneer;
mod runner;

use std::path::PathBuf;
use std::sync::Arc;

use tauri::{AppHandle, Manager, State};

struct AppState {
    jobs: Arc<runner::Jobs>,
}

/// Resolve the packaged Pioneer distribution. Called once at startup so a bad
/// install surfaces immediately rather than on the first Run click.
#[tauri::command]
fn pioneer_info(app: AppHandle) -> Result<pioneer::PioneerInfo, String> {
    let resource_dir = app.path().resource_dir().ok();
    pioneer::resolve_home(resource_dir.as_deref())
}

#[tauri::command]
fn inspect_path(path: String) -> paths::PathInfo {
    paths::inspect(&path)
}

/// What a `.poin` says about itself, for the panel under the library field.
#[tauri::command]
fn library_info(path: String) -> paths::LibraryInfo {
    paths::library_info(&path)
}

#[tauri::command]
fn read_config(path: String) -> Result<String, String> {
    paths::read_config(&path)
}

/// Show a finished run's output folder in the system file manager.
///
/// A narrow command rather than the opener plugin: all we ever need is to
/// reveal one directory we already resolved ourselves, and the plugin would
/// bring a much wider surface (arbitrary paths and URLs) for it.
///
/// The path must exist and be a directory. A run whose output has since been
/// moved or deleted therefore reports that, rather than silently doing nothing
/// or handing an arbitrary string to a shell.
#[tauri::command]
fn open_folder(path: String) -> Result<(), String> {
    let p = PathBuf::from(path.trim());
    if !p.is_dir() {
        return Err(format!("{} is no longer there.", p.display()));
    }

    #[cfg(target_os = "macos")]
    let mut cmd = std::process::Command::new("open");
    #[cfg(target_os = "windows")]
    let mut cmd = std::process::Command::new("explorer");
    #[cfg(all(unix, not(target_os = "macos")))]
    let mut cmd = std::process::Command::new("xdg-open");

    runner::hide_console(&mut cmd);
    cmd.arg(&p);
    // `explorer` returns a non-zero exit code even when it succeeds, so the
    // status is deliberately not checked -- only the spawn is.
    cmd.spawn()
        .map(|_| ())
        .map_err(|e| format!("Could not open {}: {e}", p.display()))
}

/// Number of logical cores, for the thread picker's ceiling.
#[tauri::command]
fn cpu_count() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8)
}

#[tauri::command]
fn start_job(
    app: AppHandle,
    state: State<'_, AppState>,
    job_id: String,
    command: pioneer::Command,
    invocation: runner::Invocation,
    threads: u32,
) -> Result<runner::Started, String> {
    let resource_dir = app.path().resource_dir().ok();
    let info = pioneer::resolve_home(resource_dir.as_deref())?;
    let spec = runner::Spec {
        job_id,
        command,
        home: PathBuf::from(&info.home),
        invocation,
        threads: threads.max(1),
    };
    runner::start(app.clone(), Arc::clone(&state.jobs), spec)
}

#[tauri::command]
fn cancel_job(state: State<'_, AppState>, job_id: String) -> Result<(), String> {
    state.jobs.cancel(&job_id)
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .manage(AppState { jobs: Arc::new(runner::Jobs::default()) })
        .invoke_handler(tauri::generate_handler![
            pioneer_info,
            inspect_path,
            read_config,
            library_info,
            cpu_count,
            open_folder,
            start_job,
            cancel_job,
        ])
        .run(tauri::generate_context!())
        .expect("error while running Pioneer console");
}

#[cfg(test)]
mod open_folder_tests {
    use super::open_folder;

    #[test]
    fn refuses_a_path_that_is_not_a_directory() {
        // The case a stale history entry produces: the run is still listed but
        // its output has been moved or deleted.
        let gone = std::env::temp_dir().join("pioneer_open_folder_missing");
        let _ = std::fs::remove_dir_all(&gone);
        let err = open_folder(gone.to_string_lossy().into_owned()).unwrap_err();
        assert!(err.contains("no longer there"), "unexpected message: {err}");

        // A file is not a folder either.
        let f = std::env::temp_dir().join("pioneer_open_folder_file.txt");
        std::fs::write(&f, b"x").unwrap();
        assert!(open_folder(f.to_string_lossy().into_owned()).is_err());
        let _ = std::fs::remove_file(&f);
    }
}
