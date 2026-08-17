mod history;
mod paths;
mod pioneer;
mod runner;

use std::path::PathBuf;
use std::sync::Arc;

use tauri::{AppHandle, Manager, State};

struct AppState {
    jobs: Arc<runner::Jobs>,
    /// None when the store could not be opened; the frontend then keeps using
    /// localStorage rather than losing history entirely.
    history: Option<history::Store>,
}

fn now_secs() -> i64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0)
}

fn store<'a>(state: &'a State<'_, AppState>) -> Result<&'a history::Store, String> {
    state.history.as_ref().ok_or_else(|| "history store unavailable".to_string())
}

#[tauri::command]
fn history_load(state: State<'_, AppState>) -> Result<Vec<history::Run>, String> {
    store(&state)?.load()
}

#[tauri::command]
fn history_save(state: State<'_, AppState>, run: history::Run) -> Result<(), String> {
    store(&state)?.save(&run, now_secs())
}

#[tauri::command]
fn history_delete(state: State<'_, AppState>, id: String) -> Result<(), String> {
    store(&state)?.delete(&id)
}

#[tauri::command]
fn history_next_run_no(state: State<'_, AppState>) -> Result<i64, String> {
    store(&state)?.next_run_no()
}

#[tauri::command]
fn history_needs_import(state: State<'_, AppState>) -> Result<bool, String> {
    Ok(!store(&state)?.imported()?)
}

/// Take the old localStorage rows. Called once, by the frontend, since only it
/// can read localStorage.
#[tauri::command]
fn history_import(
    state: State<'_, AppState>,
    runs: Vec<history::Run>,
    counter: i64,
) -> Result<(), String> {
    store(&state)?.import(&runs, counter, now_secs())
}

/// Resolve the packaged Pioneer distribution. Called once at startup so a bad
/// install surfaces immediately rather than on the first Run click.
#[tauri::command]
fn pioneer_info(app: AppHandle) -> Result<pioneer::PioneerInfo, String> {
    let resource_dir = app.path().resource_dir().ok();
    pioneer::resolve_home(resource_dir.as_deref())
}

/// The GUI's own version, for the sidebar footer. Compiled in rather than read
/// at runtime, so it cannot disagree with the binary the user is running.
#[tauri::command]
fn app_version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
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

/// The catalog of downloadable libraries, as JSON.
///
/// Runs `bin/DownloadSpecLib --list --json` and returns its stdout verbatim for
/// the frontend to parse. Unlike a search this is a short, quiet call whose
/// output is *data*, so it is captured rather than streamed through the job
/// runner -- the log drawer would otherwise fill with a JSON blob.
/// `async` and off-thread deliberately. Listing costs ~2.5s -- about half the
/// binary's own start-up, the rest two HTTPS round trips -- and a synchronous
/// command runs on the main thread, which froze the window for the duration and
/// meant the button's own "Loading…" state never got painted.
#[tauri::command]
async fn list_spec_libs(app: AppHandle, repo: Option<String>) -> Result<String, String> {
    let resource_dir = app.path().resource_dir().ok();
    let info = pioneer::resolve_home(resource_dir.as_deref())?;
    let resolved = pioneer::resolve_command(
        std::path::Path::new(&info.home), pioneer::Command::DownloadSpecLib)?;

    tauri::async_runtime::spawn_blocking(move || {
        let mut command = std::process::Command::new(&resolved.program);
        runner::hide_console(&mut command);
        command.args(&resolved.leading_args).arg("--list").arg("--json");
        if let Some(repo) = repo.as_deref() {
            if !repo.trim().is_empty() {
                command.arg("--repo").arg(repo.trim());
            }
        }

        let output = command.output().map_err(|e| format!(
            "could not run {}: {e}", resolved.program.display()))?;
        if !output.status.success() {
            // Pioneer writes its diagnostics to stderr; surfacing them beats a
            // bare exit code when the failure is "no network" or "repo not
            // found".
            return Err(format!("DownloadSpecLib --list failed: {}",
                               String::from_utf8_lossy(&output.stderr)));
        }
        Ok(String::from_utf8_lossy(&output.stdout).into_owned())
    })
    .await
    .map_err(|e| format!("listing task failed: {e}"))?
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
        // Must be registered first, before anything else can take effect: a
        // second launch hands its arguments to the running instance and exits.
        //
        // Without it two processes share one WebKit profile and therefore one
        // localStorage, where the whole run history is written as a single blob
        // on every change -- so whichever instance saved last silently discarded
        // the other's runs. They would also each run a queue with no knowledge
        // of the other, happily asking for 15 threads apiece on a 16-core box.
        .plugin(tauri_plugin_single_instance::init(|app, _args, _cwd| {
            if let Some(w) = app.get_webview_window("main") {
                let _ = w.unminimize();
                let _ = w.set_focus();
            }
        }))
        .plugin(tauri_plugin_dialog::init())
        .setup(|app| {
            // Opened here rather than in `manage` above because it needs the
            // AppHandle to resolve the data directory. A failure is not fatal:
            // the frontend falls back to localStorage.
            let history = match history::Store::open(&app.handle().clone()) {
                Ok(s) => Some(s),
                Err(e) => {
                    eprintln!("Pioneer console: history unavailable ({e}); using localStorage");
                    None
                }
            };
            app.manage(AppState { jobs: Arc::new(runner::Jobs::default()), history });
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![
            pioneer_info,
            app_version,
            inspect_path,
            read_config,
            library_info,
            list_spec_libs,
            cpu_count,
            open_folder,
            start_job,
            cancel_job,
            history_load,
            history_save,
            history_delete,
            history_next_run_no,
            history_needs_import,
            history_import,
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
