//! Spawning a Pioneer executable and streaming its output to the frontend.
//!
//! One job == one child process. Output is streamed line-by-line as it is
//! produced (no buffering until exit, so a multi-hour search is watchable), and
//! a job can be cancelled by killing its whole process group.

use std::collections::HashMap;
use std::io::BufReader;
use std::path::PathBuf;
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};

use tauri::{AppHandle, Emitter};

use crate::pioneer;

/// Keep spawned processes from opening a console window on Windows.
///
/// The Pioneer executables are console-subsystem programs. Started from a GUI
/// process they get a console of their own unless CREATE_NO_WINDOW is set, so
/// every run flashed up a black window alongside the app — and the cancel path
/// did it a second time for `taskkill`.
///
/// No-op off Windows, so call sites do not need their own cfg.
#[cfg(windows)]
const CREATE_NO_WINDOW: u32 = 0x0800_0000;

#[cfg(windows)]
pub fn hide_console(cmd: &mut Command) {
    use std::os::windows::process::CommandExt;
    cmd.creation_flags(CREATE_NO_WINDOW);
}

#[cfg(not(windows))]
pub fn hide_console(_cmd: &mut Command) {}

/// A line of output from a running job.
#[derive(Clone, serde::Serialize)]
pub struct LineEvent {
    pub job_id: String,
    pub line: String,
    /// "stdout" or "stderr".
    pub stream: &'static str,
    /// True when this line was terminated by a carriage return rather than a
    /// newline — i.e. the writer intends to overwrite it (a progress bar
    /// redrawing in place). The frontend replaces the previous line instead of
    /// appending, which is what a terminal does.
    pub transient: bool,
    /// How many already-committed lines this one should replace, from an
    /// `ESC[<n>A` (cursor up) seen before it. ProgressBars.jl redraws by
    /// committing a line with `\n`, moving the cursor back up, and writing over
    /// it — so without this every redraw becomes a separate line.
    pub overwrite: u32,
}

/// A job reaching a terminal state.
#[derive(Clone, serde::Serialize)]
pub struct ExitEvent {
    pub job_id: String,
    pub code: Option<i32>,
    pub success: bool,
    /// Set when the job was cancelled from the UI rather than exiting on its own.
    pub cancelled: bool,
    /// Human-readable summary for the failure banner.
    pub message: String,
}

/// Live child processes, keyed by job id.
///
/// We keep only the pid and a cancel flag here — the `Child` itself is moved
/// into the waiter thread so nothing has to share it across a lock.
#[derive(Default)]
pub struct Jobs {
    inner: Mutex<HashMap<String, JobHandle>>,
}

struct JobHandle {
    pid: u32,
    cancelled: Arc<AtomicBool>,
}

impl Jobs {
    pub fn has_active(&self) -> Result<bool, String> {
        self.inner
            .lock()
            .map(|map| !map.is_empty())
            .map_err(|e| e.to_string())
    }

    pub fn cancel(&self, job_id: &str) -> Result<(), String> {
        let map = self.inner.lock().map_err(|e| e.to_string())?;
        let handle = map
            .get(job_id)
            .ok_or_else(|| format!("job {job_id} is not running"))?;
        handle.cancelled.store(true, Ordering::SeqCst);
        kill_tree(handle.pid);
        Ok(())
    }

    fn insert(&self, job_id: String, handle: JobHandle) {
        if let Ok(mut map) = self.inner.lock() {
            map.insert(job_id, handle);
        }
    }

    fn remove(&self, job_id: &str) {
        if let Ok(mut map) = self.inner.lock() {
            map.remove(job_id);
        }
    }
}

/// How a command is configured on the command line.
///
/// The two Julia executables take a single positional argument: the path to a
/// params JSON. PioneerConverter is a separate .NET program that takes a
/// positional RAW path plus flags and has no params file at all — so the
/// invocation shape, not just the arguments, differs per command.
#[derive(serde::Deserialize)]
#[serde(tag = "kind", rename_all = "camelCase")]
pub enum Invocation {
    /// Serialize `json` to disk and pass its path.
    ParamsFile { json: String },
    /// Pass these arguments through verbatim.
    Args { args: Vec<String> },
    /// Run the same program once per argument set, in order, as one job.
    ///
    /// For PioneerConverter over a chosen list of `.raw` files: its Thermo
    /// reader takes one path and cannot read a linked file, so the list can
    /// only be one invocation each — but that is one thing the user asked for
    /// and should be one row in the queue, one log, and one outcome.
    Steps { steps: Vec<Vec<String>> },
}

/// Everything needed to start one run.
pub struct Spec {
    pub job_id: String,
    pub command: pioneer::Command,
    pub home: PathBuf,
    pub invocation: Invocation,
    pub threads: u32,
}

/// Where the params file for a job is written.
///
/// Kept after the run rather than deleted — it is the exact input Pioneer saw,
/// which is what you want when a run needs to be explained or reproduced.
fn params_path(job_id: &str, command: pioneer::Command) -> PathBuf {
    let dir = std::env::temp_dir().join("pioneer-console").join(job_id);
    let name = match command {
        pioneer::Command::SearchDia => "search.json",
        pioneer::Command::BuildSpecLib => "buildspeclib.json",
        // Argv-driven, so it never writes one; named for completeness.
        pioneer::Command::DownloadSpecLib => "download.json",
        pioneer::Command::ConvertRaw => "convert.json",
        pioneer::Command::ConvertMzml => "convertmzml.json",
    };
    dir.join(name)
}

/// What a freshly started job reports back to the UI.
#[derive(serde::Serialize)]
pub struct Started {
    pub params_path: String,
    /// The environment we actually set, so the log shows what ran rather than
    /// the frontend guessing at the same formula.
    pub env_summary: String,
}

/// Spawn the job and stream its output as events.
///
/// A job is usually one child process, but not always: `Invocation::Steps` runs
/// several in order under a single job id, streaming into one log and reporting
/// one outcome. That exists because PioneerConverter's Thermo reader cannot be
/// handed a list and cannot read a linked file, so converting a chosen set of
/// `.raw` means one invocation per file — while the queue should still show the
/// one thing the user asked for.
pub fn start(app: AppHandle, jobs: Arc<Jobs>, spec: Spec) -> Result<Started, String> {
    let resolved = pioneer::resolve_command(&spec.home, spec.command)?;

    // `params_path` doubles as the record of what was run: for a params-file
    // command it is the JSON; for an argv command there is no file, so the
    // reported path is empty and the log carries the command line instead.
    let mut params_display = String::new();
    let steps: Vec<Vec<String>> = match &spec.invocation {
        Invocation::ParamsFile { json } => {
            let params = params_path(&spec.job_id, spec.command);
            if let Some(parent) = params.parent() {
                std::fs::create_dir_all(parent)
                    .map_err(|e| format!("could not create {}: {e}", parent.display()))?;
            }
            std::fs::write(&params, json)
                .map_err(|e| format!("could not write {}: {e}", params.display()))?;
            params_display = params.display().to_string();
            vec![vec![params.display().to_string()]]
        }
        Invocation::Args { args } => vec![args.clone()],
        Invocation::Steps { steps } => steps.clone(),
    };
    if steps.is_empty() {
        return Err("nothing to run".to_string());
    }

    // Bypassing the `pioneer` wrapper means we own these. Without
    // JULIA_NUM_THREADS a direct bin/SearchDIA call defaults to a single thread.
    let mut env_summary = String::new();
    let mut envs: Vec<(String, String)> = Vec::new();
    if spec.command.is_julia() && !resolved.via_wrapper {
        // `JULIA_NUM_GC_THREADS` is the variable Julia actually reads — it is
        // the env form of `--gcthreads=N,M` (N marking, M concurrent sweeper).
        // `JULIA_GC_THREADS` is not a Julia variable and is silently ignored,
        // which leaves GC threads defaulting to the full worker count.
        // Matches the `pioneer` wrapper's own formula, `(threads + 1) / 2`,
        // i.e. ceil rather than floor — so the GUI and the shell script agree.
        let gc = format!("{},1", ((spec.threads + 1) / 2).max(1));
        envs.push(("JULIA_NUM_THREADS".into(), spec.threads.to_string()));
        envs.push(("JULIA_NUM_GC_THREADS".into(), gc.clone()));
        env_summary = format!(
            "JULIA_NUM_THREADS={} JULIA_NUM_GC_THREADS={}",
            spec.threads, gc
        );
    }

    // The first step is spawned here rather than in the driver thread so that a
    // distribution that cannot be executed at all still fails out of `start`,
    // where the frontend surfaces it as an error on the Run button.
    let multi = steps.len() > 1;
    if multi {
        announce(&app, &spec.job_id, &resolved, &steps[0], 1, steps.len());
    }
    let first = spawn_step(&resolved, &steps[0], &envs)?;

    let cancelled = Arc::new(AtomicBool::new(false));
    jobs.insert(
        spec.job_id.clone(),
        JobHandle { pid: first.id(), cancelled: Arc::clone(&cancelled) },
    );

    // Driver thread: owns each Child in turn, and reports one terminal state for
    // the whole sequence.
    let job_id = spec.job_id.clone();
    let jobs_for_wait = Arc::clone(&jobs);
    std::thread::spawn(move || {
        let mut child = first;
        let mut index = 0usize;
        let event = loop {
            // Drain both pipes to completion before the next step starts, so a
            // step's output cannot interleave with the tail of the one before.
            let out = child.stdout.take().map(|r| pump(app.clone(), job_id.clone(), r, "stdout"));
            let err = child.stderr.take().map(|r| pump(app.clone(), job_id.clone(), r, "stderr"));
            let status = child.wait();
            if let Some(h) = out {
                let _ = h.join();
            }
            if let Some(h) = err {
                let _ = h.join();
            }
            let was_cancelled = cancelled.load(Ordering::SeqCst);

            let st = match status {
                Ok(st) => st,
                Err(e) => {
                    break ExitEvent {
                        job_id: job_id.clone(),
                        code: None,
                        success: false,
                        cancelled: was_cancelled,
                        message: format!("Could not wait for Pioneer: {e}"),
                    }
                }
            };

            index += 1;
            let last = index >= steps.len();
            // A failed step ends the job: the steps of a sequence are the same
            // command over different files, and continuing after one has failed
            // would bury the failure under whatever came next.
            if was_cancelled || !st.success() || last {
                let success = st.success() && !was_cancelled;
                let code = st.code();
                let message = if was_cancelled {
                    "Cancelled by user.".to_string()
                } else if success {
                    String::new()
                } else {
                    let where_ = if multi {
                        format!(" on step {index} of {}", steps.len())
                    } else {
                        String::new()
                    };
                    match code {
                        Some(c) => format!("Pioneer exited with code {c}{where_}."),
                        None => format!("Pioneer was terminated by a signal{where_}."),
                    }
                };
                break ExitEvent {
                    job_id: job_id.clone(),
                    code,
                    success,
                    cancelled: was_cancelled,
                    message,
                };
            }

            announce(&app, &job_id, &resolved, &steps[index], index + 1, steps.len());
            match spawn_step(&resolved, &steps[index], &envs) {
                Ok(next) => {
                    // Replaces the previous entry, so a cancel now kills the
                    // step that is actually running.
                    jobs_for_wait.insert(
                        job_id.clone(),
                        JobHandle { pid: next.id(), cancelled: Arc::clone(&cancelled) },
                    );
                    child = next;
                }
                Err(e) => {
                    break ExitEvent {
                        job_id: job_id.clone(),
                        code: None,
                        success: false,
                        cancelled: false,
                        message: e,
                    }
                }
            }
        };
        jobs_for_wait.remove(&job_id);
        let _ = app.emit("job-exit", event);
    });

    Ok(Started { params_path: params_display, env_summary })
}

/// Build and spawn one step of a job.
fn spawn_step(
    resolved: &pioneer::Resolved,
    args: &[String],
    envs: &[(String, String)],
) -> Result<Child, String> {
    let mut cmd = Command::new(&resolved.program);
    hide_console(&mut cmd);
    for arg in &resolved.leading_args {
        cmd.arg(arg);
    }
    for a in args {
        cmd.arg(a);
    }
    cmd.stdout(Stdio::piped()).stderr(Stdio::piped()).stdin(Stdio::null());
    for (k, v) in envs {
        cmd.env(k, v);
    }
    // Put the child in its own process group so cancel can take down Julia's
    // worker processes too, not just the parent.
    #[cfg(unix)]
    {
        use std::os::unix::process::CommandExt;
        cmd.process_group(0);
    }
    cmd.spawn()
        .map_err(|e| format!("could not start {}: {e}", resolved.program.display()))
}

/// Write a header into the log naming the step about to run.
///
/// Only for multi-step jobs, where the log would otherwise be several files'
/// output run together with nothing saying where one ends and the next begins.
fn announce(
    app: &AppHandle,
    job_id: &str,
    resolved: &pioneer::Resolved,
    args: &[String],
    n: usize,
    total: usize,
) {
    let program = resolved
        .program
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| resolved.program.display().to_string());
    let quoted: Vec<String> = args
        .iter()
        .map(|a| if a.contains(' ') { format!("\"{a}\"") } else { a.clone() })
        .collect();
    let _ = app.emit(
        "job-line",
        LineEvent {
            job_id: job_id.to_string(),
            line: format!("[{n}/{total}] {program} {}", quoted.join(" ")),
            stream: "stdout",
            transient: false,
            overwrite: 0,
        },
    );
}

/// Forward one pipe to the frontend, honouring terminal carriage-return
/// semantics.
///
/// Julia's progress bars (ProgressBars.jl) redraw by emitting `\r` and
/// rewriting the line, and only send `\n` when the bar is finished. Splitting
/// on `\n` alone would either buffer a whole stage into one line or leave stray
/// CRs in the text, which the log pane renders as extra line breaks. So we
/// treat `\r` as "this line is about to be overwritten" and tell the frontend
/// to replace rather than append.
fn pump<R: std::io::Read + Send + 'static>(
    app: AppHandle,
    job_id: String,
    reader: R,
    stream: &'static str,
) -> std::thread::JoinHandle<()> {
    std::thread::spawn(move || {
        let mut buf = BufReader::new(reader);
        let mut splitter = LineSplitter::default();
        let mut byte = [0u8; 1];
        loop {
            match std::io::Read::read(&mut buf, &mut byte) {
                Ok(0) => break,
                Ok(_) => {}
                Err(_) => break,
            }
            if let Some(seg) = splitter.push(byte[0]) {
                let _ = app.emit(
                    "job-line",
                    LineEvent {
                        job_id: job_id.clone(),
                        line: seg.text,
                        stream,
                        transient: seg.transient,
                        overwrite: seg.overwrite,
                    },
                );
            }
        }
        // Whatever was mid-line when the pipe closed.
        if let Some(seg) = splitter.finish() {
            let _ = app.emit(
                "job-line",
                LineEvent {
                    job_id: job_id.clone(),
                    line: seg.text,
                    stream,
                    transient: seg.transient,
                    overwrite: seg.overwrite,
                },
            );
        }
    })
}

/// One line ready to show, with the terminal control that applied to it.
pub struct Segment {
    pub text: String,
    pub transient: bool,
    pub overwrite: u32,
}

/// Where the little ANSI parser is in an `ESC [ params final` sequence.
enum Ansi {
    Normal,
    Esc,
    Csi,
}

impl Default for Ansi {
    fn default() -> Self {
        Ansi::Normal
    }
}

/// Splits a byte stream into displayable lines, interpreting the small amount
/// of terminal control Julia's progress bars use: carriage returns and
/// `ESC[<n>A` (cursor up). Everything else in a CSI sequence is discarded.
#[derive(Default)]
pub struct LineSplitter {
    cur: Vec<u8>,
    /// A CR is only an overwrite if it is not part of a CRLF pair, so the
    /// decision waits for the following byte.
    pending_cr: bool,
    /// Lines the next segment should replace, from cursor-up codes.
    pending_up: u32,
    ansi: Ansi,
    csi_params: Vec<u8>,
}

impl LineSplitter {
    fn take(&mut self, transient: bool) -> Option<Segment> {
        // from_utf8_lossy so a partially-written multibyte glyph in a redrawn
        // bar degrades instead of killing the stream.
        let text = String::from_utf8_lossy(&self.cur).to_string();
        self.cur.clear();
        let overwrite = std::mem::take(&mut self.pending_up);
        Some(Segment { text, transient, overwrite })
    }

    pub fn push(&mut self, b: u8) -> Option<Segment> {
        match self.ansi {
            Ansi::Esc => {
                self.ansi = if b == b'[' {
                    self.csi_params.clear();
                    Ansi::Csi
                } else {
                    Ansi::Normal
                };
                return None;
            }
            Ansi::Csi => {
                // The final byte of a CSI sequence is in 0x40..=0x7E.
                if (0x40..=0x7e).contains(&b) {
                    if b == b'A' {
                        let n = std::str::from_utf8(&self.csi_params)
                            .ok()
                            .and_then(|s| s.trim().parse::<u32>().ok())
                            .unwrap_or(1);
                        self.pending_up += n.max(1);
                    }
                    // Every other CSI code (colour, erase, …) is dropped.
                    self.ansi = Ansi::Normal;
                } else {
                    self.csi_params.push(b);
                }
                return None;
            }
            Ansi::Normal => {
                if b == 0x1b {
                    self.ansi = Ansi::Esc;
                    return None;
                }
            }
        }

        if self.pending_cr {
            self.pending_cr = false;
            if b == b'\n' {
                // CRLF — an ordinary line ending.
                return self.take(false);
            }
            // A bare CR: commit what we have as overwritable. Skip empties so a
            // leading "\r" does not produce a blank line.
            if !self.cur.is_empty() {
                let seg = self.take(true);
                // The byte that ended the CR still belongs to the next line.
                if b != b'\r' {
                    self.cur.push(b);
                } else {
                    self.pending_cr = true;
                }
                return seg;
            }
        }

        match b {
            b'\n' => self.take(false),
            b'\r' => {
                self.pending_cr = true;
                None
            }
            _ => {
                self.cur.push(b);
                None
            }
        }
    }

    pub fn finish(&mut self) -> Option<Segment> {
        if self.cur.is_empty() {
            return None;
        }
        self.take(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Feed bytes through the splitter and apply the frontend's replace rule,
    /// returning what the log pane would end up showing.
    fn render(input: &[u8]) -> Vec<String> {
        let mut splitter = LineSplitter::default();
        let mut segs: Vec<Segment> = Vec::new();
        for &b in input {
            if let Some(s) = splitter.push(b) {
                segs.push(s);
            }
        }
        if let Some(s) = splitter.finish() {
            segs.push(s);
        }

        let mut lines: Vec<(String, bool)> = Vec::new();
        for s in segs {
            let mut drop = s.overwrite;
            while drop > 0 && !lines.is_empty() {
                lines.pop();
                drop -= 1;
            }
            if s.overwrite == 0 {
                if let Some(last) = lines.last() {
                    if last.1 {
                        lines.pop();
                    }
                }
            }
            lines.push((s.text, s.transient));
        }
        lines.into_iter().map(|(t, _)| t).collect()
    }

    #[test]
    fn plain_lines_pass_through() {
        assert_eq!(render(b"one\ntwo\nthree\n"), vec!["one", "two", "three"]);
    }

    #[test]
    fn crlf_is_an_ordinary_line_ending() {
        assert_eq!(render(b"one\r\ntwo\r\n"), vec!["one", "two"]);
    }

    /// The shape ProgressBars.jl actually emits: commit with `\n`, then move
    /// the cursor back up and redraw. All redraws must collapse to one line.
    #[test]
    fn progress_bar_redraw_collapses_to_one_line() {
        let input = b"\r0% bar\n\x1b[1A\r50% bar\n\x1b[1A\r100% bar\n";
        assert_eq!(render(input), vec!["100% bar"]);
    }

    /// A bar sequence followed by ordinary output must leave the finished bar
    /// in place — the cursor-up applies to the bar, not to what came after.
    #[test]
    fn output_after_a_bar_is_kept() {
        let input = b"\r0% bar\n\x1b[1A\r100% bar\ndone\n";
        assert_eq!(render(input), vec!["100% bar", "done"]);
    }

    /// A writer that redraws with bare carriage returns and never commits.
    #[test]
    fn bare_cr_redraw_collapses() {
        assert_eq!(render(b"10%\r20%\r30%\n"), vec!["30%"]);
    }

    #[test]
    fn non_cursor_ansi_codes_are_stripped() {
        assert_eq!(render(b"\x1b[31mred\x1b[0m\n"), vec!["red"]);
    }

    #[test]
    fn multi_line_cursor_up_removes_several() {
        let input = b"a\nb\nc\n\x1b[2Ax\n";
        assert_eq!(render(input), vec!["a", "x"]);
    }
}

/// Kill a process group (unix) or process tree (windows).
fn kill_tree(pid: u32) {
    #[cfg(unix)]
    {
        // Negative pid targets the process group established via process_group(0).
        let _ = Command::new("/bin/kill")
            .args(["-TERM", &format!("-{pid}")])
            .status();
        let pid_copy = pid;
        std::thread::spawn(move || {
            std::thread::sleep(std::time::Duration::from_secs(5));
            let _ = Command::new("/bin/kill")
                .args(["-KILL", &format!("-{pid_copy}")])
                .status();
        });
    }
    #[cfg(windows)]
    {
        let mut kill = Command::new("taskkill");
        hide_console(&mut kill);
        let _ = kill.args(["/PID", &pid.to_string(), "/T", "/F"]).status();
    }
}
