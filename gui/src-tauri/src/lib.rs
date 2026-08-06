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

#[tauri::command]
fn read_config(path: String) -> Result<String, String> {
    paths::read_config(&path)
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
            cpu_count,
            start_job,
            cancel_job,
        ])
        .run(tauri::generate_context!())
        .expect("error while running Pioneer console");
}
