//! Run history, kept in SQLite under the app's data directory.
//!
//! Replaces a `localStorage` array that was rewritten whole on every change.
//! That had three problems this does not: it was capped to keep the blob small,
//! a quota failure could lose the lot, and two instances writing the same blob
//! silently discarded each other's runs.
//!
//! The database lives beside the app's other data rather than inside the
//! install directory, so reinstalling or upgrading does not take the history
//! with it:
//!
//!   macOS    ~/Library/Application Support/org.pioneer.console/history.db
//!   Windows  %APPDATA%\org.pioneer.console\history.db
//!   Linux    ~/.local/share/org.pioneer.console/history.db
//!
//! WAL mode, so a reader never blocks a writer. The app is single-instance, but
//! the database is a file on disk and something else may well read it one day.

use std::path::PathBuf;
use std::sync::Mutex;

use rusqlite::{params, Connection};

/// One finished run, as the frontend stores and restores it.
///
/// `snapshot` is the form state as JSON. It is opaque here on purpose: the
/// three commands have quite different parameter sets, and teaching this layer
/// about them would mean changing it every time a form field is added.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Run {
    pub id: String,
    pub run_no: i64,
    pub cmd: String,
    pub title: String,
    pub target: String,
    pub threads: i64,
    pub status: String,
    pub snapshot: String,
    /// Unix seconds. Set when the row is written.
    #[serde(default)]
    pub finished_at: i64,
}

pub struct Store {
    conn: Mutex<Connection>,
}

fn data_dir(app: &tauri::AppHandle) -> Result<PathBuf, String> {
    use tauri::Manager;
    app.path()
        .app_data_dir()
        .map_err(|e| format!("No app data directory: {e}"))
}

impl Store {
    pub fn open(app: &tauri::AppHandle) -> Result<Self, String> {
        let dir = data_dir(app)?;
        std::fs::create_dir_all(&dir).map_err(|e| format!("Could not create {dir:?}: {e}"))?;
        let conn = Connection::open(dir.join("history.db"))
            .map_err(|e| format!("Could not open history.db: {e}"))?;
        // WAL survives a crash mid-write and lets a reader run during one.
        conn.pragma_update(None, "journal_mode", "WAL")
            .map_err(|e| format!("Could not enable WAL: {e}"))?;
        conn.execute_batch(
            "CREATE TABLE IF NOT EXISTS runs (
                 id          TEXT PRIMARY KEY,
                 run_no      INTEGER NOT NULL,
                 cmd         TEXT NOT NULL,
                 title       TEXT NOT NULL,
                 target      TEXT NOT NULL,
                 threads     INTEGER NOT NULL,
                 status      TEXT NOT NULL,
                 snapshot    TEXT NOT NULL,
                 finished_at INTEGER NOT NULL
             );
             CREATE INDEX IF NOT EXISTS runs_by_number ON runs(run_no DESC);
             CREATE TABLE IF NOT EXISTS meta (
                 key   TEXT PRIMARY KEY,
                 value TEXT NOT NULL
             );",
        )
        .map_err(|e| format!("Could not create the schema: {e}"))?;
        Ok(Store { conn: Mutex::new(conn) })
    }

    fn lock(&self) -> Result<std::sync::MutexGuard<'_, Connection>, String> {
        self.conn.lock().map_err(|_| "history lock poisoned".to_string())
    }

    /// Oldest first, matching the order the sidebar renders.
    pub fn load(&self) -> Result<Vec<Run>, String> {
        let conn = self.lock()?;
        let mut stmt = conn
            .prepare(
                "SELECT id, run_no, cmd, title, target, threads, status, snapshot, finished_at
                 FROM runs ORDER BY run_no ASC",
            )
            .map_err(|e| e.to_string())?;
        let rows = stmt
            .query_map([], |r| {
                Ok(Run {
                    id: r.get(0)?,
                    run_no: r.get(1)?,
                    cmd: r.get(2)?,
                    title: r.get(3)?,
                    target: r.get(4)?,
                    threads: r.get(5)?,
                    status: r.get(6)?,
                    snapshot: r.get(7)?,
                    finished_at: r.get(8)?,
                })
            })
            .map_err(|e| e.to_string())?;
        rows.collect::<Result<Vec<_>, _>>().map_err(|e| e.to_string())
    }

    /// Insert or update one run. Called as each finishes, not for the whole
    /// list — the point of moving off localStorage is not rewriting everything.
    pub fn save(&self, run: &Run, now: i64) -> Result<(), String> {
        let conn = self.lock()?;
        conn.execute(
            "INSERT INTO runs (id, run_no, cmd, title, target, threads, status, snapshot, finished_at)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)
             ON CONFLICT(id) DO UPDATE SET
                 status = excluded.status,
                 target = excluded.target,
                 snapshot = excluded.snapshot",
            params![
                run.id, run.run_no, run.cmd, run.title, run.target, run.threads, run.status,
                run.snapshot,
                if run.finished_at > 0 { run.finished_at } else { now },
            ],
        )
        .map_err(|e| e.to_string())?;
        Ok(())
    }

    pub fn delete(&self, id: &str) -> Result<(), String> {
        let conn = self.lock()?;
        conn.execute("DELETE FROM runs WHERE id = ?1", params![id])
            .map_err(|e| e.to_string())?;
        Ok(())
    }

    /// Hand out the next history number.
    ///
    /// Kept in the database rather than derived from `MAX(run_no)`, which would
    /// rewind as soon as the newest run was deleted — the numbers are an
    /// identity, not a position.
    pub fn next_run_no(&self) -> Result<i64, String> {
        let conn = self.lock()?;
        conn.execute(
            "INSERT INTO meta (key, value) VALUES ('run_counter', '1')
             ON CONFLICT(key) DO UPDATE SET value = CAST(CAST(value AS INTEGER) + 1 AS TEXT)",
            [],
        )
        .map_err(|e| e.to_string())?;
        let v: String = conn
            .query_row("SELECT value FROM meta WHERE key = 'run_counter'", [], |r| r.get(0))
            .map_err(|e| e.to_string())?;
        v.parse::<i64>().map_err(|e| e.to_string())
    }

    /// True once the localStorage import has been done, so it happens once.
    pub fn imported(&self) -> Result<bool, String> {
        let conn = self.lock()?;
        let n: i64 = conn
            .query_row("SELECT COUNT(*) FROM meta WHERE key = 'imported_v1'", [], |r| r.get(0))
            .map_err(|e| e.to_string())?;
        Ok(n > 0)
    }

    /// Import the old localStorage rows, and remember that we did.
    ///
    /// The counter is carried across too: without it the first run after an
    /// upgrade would be numbered 1 again, colliding with an existing entry.
    pub fn import(&self, runs: &[Run], counter: i64, now: i64) -> Result<(), String> {
        for r in runs {
            self.save(r, now)?;
        }
        let conn = self.lock()?;
        let highest = runs.iter().map(|r| r.run_no).max().unwrap_or(0).max(counter);
        conn.execute(
            "INSERT OR REPLACE INTO meta (key, value) VALUES ('run_counter', ?1)",
            params![highest.to_string()],
        )
        .map_err(|e| e.to_string())?;
        conn.execute(
            "INSERT OR REPLACE INTO meta (key, value) VALUES ('imported_v1', '1')",
            [],
        )
        .map_err(|e| e.to_string())?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn store() -> Store {
        let conn = Connection::open_in_memory().unwrap();
        conn.execute_batch(
            "CREATE TABLE runs (id TEXT PRIMARY KEY, run_no INTEGER NOT NULL, cmd TEXT NOT NULL,
                 title TEXT NOT NULL, target TEXT NOT NULL, threads INTEGER NOT NULL,
                 status TEXT NOT NULL, snapshot TEXT NOT NULL, finished_at INTEGER NOT NULL);
             CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);",
        )
        .unwrap();
        Store { conn: Mutex::new(conn) }
    }

    fn run(id: &str, no: i64) -> Run {
        Run {
            id: id.into(),
            run_no: no,
            cmd: "searchdia".into(),
            title: "crisp-bison".into(),
            target: "/tmp/out".into(),
            threads: 15,
            status: "done".into(),
            snapshot: "{}".into(),
            finished_at: 0,
        }
    }

    #[test]
    fn saves_and_reloads_in_run_order() {
        let s = store();
        s.save(&run("b", 2), 100).unwrap();
        s.save(&run("a", 1), 100).unwrap();
        let got = s.load().unwrap();
        assert_eq!(got.iter().map(|r| r.run_no).collect::<Vec<_>>(), vec![1, 2]);
    }

    #[test]
    fn saving_the_same_id_updates_rather_than_duplicating() {
        let s = store();
        s.save(&run("a", 1), 100).unwrap();
        let mut later = run("a", 1);
        later.status = "failed".into();
        s.save(&later, 200).unwrap();
        let got = s.load().unwrap();
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].status, "failed");
        // The original timestamp stands: the run finished once.
        assert_eq!(got[0].finished_at, 100);
    }

    #[test]
    fn numbers_keep_climbing_when_the_newest_run_is_deleted() {
        // The behaviour MAX(run_no) would get wrong, and the reason the counter
        // is stored rather than derived.
        let s = store();
        assert_eq!(s.next_run_no().unwrap(), 1);
        assert_eq!(s.next_run_no().unwrap(), 2);
        s.save(&run("b", 2), 100).unwrap();
        s.delete("b").unwrap();
        assert_eq!(s.next_run_no().unwrap(), 3, "deleting must not rewind the counter");
    }

    #[test]
    fn import_runs_once_and_carries_the_counter_forward() {
        let s = store();
        assert!(!s.imported().unwrap());
        s.import(&[run("a", 7), run("b", 9)], 9, 100).unwrap();
        assert!(s.imported().unwrap());
        assert_eq!(s.load().unwrap().len(), 2);
        // Not 1: a fresh counter would collide with the runs just imported.
        assert_eq!(s.next_run_no().unwrap(), 10);
    }

    #[test]
    fn there_is_no_cap() {
        let s = store();
        for i in 1..=500 {
            s.save(&run(&format!("r{i}"), i), 100).unwrap();
        }
        assert_eq!(s.load().unwrap().len(), 500);
    }
}
