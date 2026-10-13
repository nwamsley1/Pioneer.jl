//! What kind of drive a path is on, so a search can warn before reading MS data
//! from a mechanical or optical drive.
//!
//! SearchDIA reads each MS data file out of order; on a spinning disk every such
//! read is a seek (one lab PC took 32 h for a search that should take 5-6 h).
//! The same check runs in SearchDIA itself (src/utils/storageType.jl), which
//! prints a warning; here it lets the GUI ask before starting. Never fails:
//! anything that cannot be determined is "unknown", which is never warned about.

use std::io::Read;
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

#[derive(Clone, serde::Serialize)]
pub struct StorageInfo {
    pub path: String,
    /// "ssd", "hdd", "optical", "network" or "unknown".
    pub kind: &'static str,
}

const PROBE_TIMEOUT: Duration = Duration::from_secs(10);

/// The kind of storage each path is on, probing each volume once.
pub fn kinds(paths: &[String]) -> Vec<StorageInfo> {
    let mut seen: Vec<(String, &'static str)> = Vec::new();
    let mut out = Vec::new();
    for p in paths {
        // A folder that does not exist yet (an output folder before the run) is judged by
        // the nearest existing folder above it: that is the drive it will be created on.
        let Some(real) = Path::new(p).ancestors().find_map(|a| std::fs::canonicalize(a).ok()) else {
            out.push(StorageInfo { path: p.clone(), kind: "unknown" });
            continue;
        };
        let key = volume_key(&real);
        let kind = match seen.iter().find(|(k, _)| *k == key) {
            Some((_, kind)) => *kind,
            None => {
                let kind = kind_of(&real);
                seen.push((key, kind));
                kind
            }
        };
        out.push(StorageInfo { path: p.clone(), kind });
    }
    out
}

/// Run a command and return its stdout, or None if it fails or exceeds the timeout.
fn read_with_timeout(mut cmd: Command) -> Option<String> {
    crate::runner::hide_console(&mut cmd);
    let mut child = cmd.stdout(Stdio::piped()).stderr(Stdio::null()).stdin(Stdio::null()).spawn().ok()?;
    let start = Instant::now();
    loop {
        match child.try_wait() {
            Ok(Some(status)) => {
                let mut s = String::new();
                child.stdout.take()?.read_to_string(&mut s).ok()?;
                return status.success().then_some(s);
            }
            Ok(None) if start.elapsed() < PROBE_TIMEOUT => std::thread::sleep(Duration::from_millis(20)),
            _ => {
                let _ = child.kill();
                return None;
            }
        }
    }
}

#[cfg(target_os = "windows")]
fn volume_key(p: &Path) -> String {
    let s = p.to_string_lossy().replace("\\\\?\\UNC\\", "\\\\").replace("\\\\?\\", "");
    if s.len() >= 2 && s.as_bytes()[1] == b':' {
        return s[..2].to_uppercase();
    }
    s.split('\\').filter(|x| !x.is_empty()).take(2).collect::<Vec<_>>().join("\\").to_uppercase()
}

#[cfg(unix)]
fn volume_key(p: &Path) -> String {
    use std::os::unix::fs::MetadataExt;
    std::fs::metadata(p).map(|m| m.dev().to_string()).unwrap_or_default()
}

// ── macOS ──────────────────────────────────────────────────────────────────

#[cfg(target_os = "macos")]
fn kind_of(p: &Path) -> &'static str {
    let mut df = Command::new("df");
    df.arg("-P").arg(p);
    let Some(out) = read_with_timeout(df) else { return "unknown" };
    let Some(last) = out.trim().lines().last() else { return "unknown" };
    let source = last.split_whitespace().next().unwrap_or("");
    if source.starts_with("//") || (!source.starts_with('/') && source.contains(":/")) {
        return "network";
    }
    if !source.starts_with("/dev/") {
        return "unknown";
    }
    let mut du = Command::new("diskutil");
    du.args(["info", "-plist", source]);
    read_with_timeout(du).map(|plist| macos_plist_kind(&plist)).unwrap_or("unknown")
}

/// Kind from `diskutil info -plist`: optical media first, then SolidState.
#[allow(dead_code)]
fn macos_plist_kind(plist: &str) -> &'static str {
    if plist.contains("<key>OpticalMediaType</key>") || plist.contains("<key>OpticalDeviceType</key>") {
        return "optical";
    }
    let Some(i) = plist.find("<key>SolidState</key>") else { return "unknown" };
    let rest = plist[i + "<key>SolidState</key>".len()..].trim_start();
    if rest.starts_with("<true") {
        "ssd"
    } else if rest.starts_with("<false") {
        "hdd"
    } else {
        "unknown"
    }
}

// ── Linux ──────────────────────────────────────────────────────────────────

#[allow(dead_code)]
const LINUX_NETWORK_FS: &[&str] = &[
    "nfs", "nfs4", "cifs", "smb3", "smbfs", "sshfs", "fuse.sshfs", "9p", "afs", "glusterfs",
    "fuse.glusterfs", "ceph", "lustre", "gpfs", "beegfs", "fuse.s3fs",
];

/// Undo mountinfo's octal escapes (`\040` for a space).
#[allow(dead_code)]
fn unescape_mountinfo(s: &str) -> String {
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'\\' && i + 3 < b.len() && b[i + 1..i + 4].iter().all(|c| (b'0'..=b'7').contains(c)) {
            out.push(((b[i + 1] - b'0') << 6) | ((b[i + 2] - b'0') << 3) | (b[i + 3] - b'0'));
            i += 4;
        } else {
            out.push(b[i]);
            i += 1;
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// (fstype, source) of the mount containing `path`: the longest mount point that prefixes it.
#[allow(dead_code)]
fn linux_mount_for(path: &str, mountinfo: &str) -> Option<(String, String)> {
    let mut best: Option<(usize, String, String)> = None;
    for line in mountinfo.lines() {
        let Some((pre, post)) = line.split_once(" - ") else { continue };
        let pre: Vec<&str> = pre.split_whitespace().collect();
        let post: Vec<&str> = post.split_whitespace().collect();
        if pre.len() < 5 || post.len() < 2 {
            continue;
        }
        let mp = unescape_mountinfo(pre[4]);
        let inside = path == mp || mp == "/" || path.starts_with(&format!("{}/", mp.trim_end_matches('/')));
        if inside && best.as_ref().map_or(true, |(len, _, _)| mp.len() > *len) {
            best = Some((mp.len(), post[0].to_string(), unescape_mountinfo(post[1])));
        }
    }
    best.map(|(_, fs, src)| (fs, src))
}

#[cfg(target_os = "linux")]
fn linux_block_kind(sysdev: &Path, depth: usize) -> &'static str {
    if depth > 8 {
        return "unknown";
    }
    let Ok(dev) = std::fs::canonicalize(sysdev) else { return "unknown" };
    let disk = if dev.join("partition").is_file() { dev.parent().unwrap_or(&dev).to_path_buf() } else { dev };
    let name = disk.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    if name.starts_with("sr") {
        return "optical";
    }
    if std::fs::read_to_string(disk.join("device/type")).map(|t| t.trim() == "5").unwrap_or(false) {
        return "optical";
    }
    if let Ok(entries) = std::fs::read_dir(disk.join("slaves")) {
        let kinds: Vec<&str> = entries.flatten().map(|e| linux_block_kind(&e.path(), depth + 1)).collect();
        if !kinds.is_empty() {
            if kinds.contains(&"hdd") {
                return "hdd";
            }
            if kinds.contains(&"optical") {
                return "optical";
            }
            return if kinds.iter().all(|k| *k == "ssd") { "ssd" } else { "unknown" };
        }
    }
    if name.starts_with("loop") || name.starts_with("zram") || name.starts_with("ram") {
        return "unknown";
    }
    match std::fs::read_to_string(disk.join("queue/rotational")).map(|s| s.trim().to_string()) {
        Ok(s) if s == "1" => "hdd",
        Ok(s) if s == "0" => "ssd",
        _ => "unknown",
    }
}

#[cfg(target_os = "linux")]
fn kind_of(p: &Path) -> &'static str {
    use std::os::unix::fs::MetadataExt;
    let path = p.to_string_lossy();
    let mount = std::fs::read_to_string("/proc/self/mountinfo").ok().and_then(|mi| linux_mount_for(&path, &mi));
    if let Some((fs, src)) = &mount {
        if LINUX_NETWORK_FS.contains(&fs.as_str()) || src.starts_with("//") {
            return "network";
        }
    }
    let major_minor = |dev: u64| (((dev >> 8) & 0xfff) | ((dev >> 32) & !0xfff), (dev & 0xff) | ((dev >> 12) & !0xff));
    let Ok(meta) = std::fs::metadata(p) else { return "unknown" };
    let (mut major, mut minor) = major_minor(meta.dev());
    if major == 0 {
        if let Some((_, src)) = &mount {
            if let Ok(m) = std::fs::metadata(src) {
                (major, minor) = major_minor(m.rdev());
            }
        }
    }
    if major == 0 {
        return "unknown";
    }
    linux_block_kind(Path::new(&format!("/sys/dev/block/{major}:{minor}")), 0)
}

// ── Windows ────────────────────────────────────────────────────────────────

#[cfg(target_os = "windows")]
fn kind_of(p: &Path) -> &'static str {
    let s = p.to_string_lossy().replace("\\\\?\\UNC\\", "\\\\").replace("\\\\?\\", "");
    if s.starts_with("\\\\") {
        return "network";
    }
    let Some(letter) = s.chars().next().filter(|c| c.is_ascii_alphabetic()) else { return "unknown" };
    let letter = letter.to_ascii_uppercase();
    let script = format!(
        "$ErrorActionPreference = 'Stop'\n\
         $ld = Get-CimInstance Win32_LogicalDisk -Filter \"DeviceID='{letter}:'\"\n\
         if ($ld.DriveType -eq 4) {{ 'NETWORK'; exit }}\n\
         if ($ld.DriveType -eq 5) {{ 'OPTICAL'; exit }}\n\
         $disk = Get-Partition -DriveLetter {letter} | Get-Disk\n\
         $pd = Get-PhysicalDisk | Where-Object {{ $_.DeviceId -eq [string]$disk.Number }}\n\
         if ($pd) {{ [string]$pd.MediaType }} else {{ 'Unspecified' }}"
    );
    let mut ps = Command::new("powershell");
    ps.args(["-NoProfile", "-NonInteractive", "-ExecutionPolicy", "Bypass", "-Command", &script]);
    read_with_timeout(ps).map(|w| windows_probe_kind(&w)).unwrap_or("unknown")
}

/// Kind from the word the PowerShell probe prints.
#[allow(dead_code)]
fn windows_probe_kind(word: &str) -> &'static str {
    match word.trim().to_uppercase().as_str() {
        "NETWORK" => "network",
        "OPTICAL" => "optical",
        "HDD" => "hdd",
        "SSD" | "SCM" => "ssd",
        _ => "unknown",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_diskutil_plists() {
        assert_eq!(macos_plist_kind("<key>SolidState</key>\n\t<true/>"), "ssd");
        assert_eq!(macos_plist_kind("<key>SolidState</key>\n\t<false/>"), "hdd");
        assert_eq!(macos_plist_kind("<key>OpticalMediaType</key><string>DVD-ROM</string>"), "optical");
        assert_eq!(macos_plist_kind("<key>Internal</key><true/>"), "unknown");
    }

    #[test]
    fn parses_windows_probe_words() {
        assert_eq!(windows_probe_kind("HDD\r\n"), "hdd");
        assert_eq!(windows_probe_kind("SSD"), "ssd");
        assert_eq!(windows_probe_kind("OPTICAL"), "optical");
        assert_eq!(windows_probe_kind("NETWORK"), "network");
        assert_eq!(windows_probe_kind("Unspecified"), "unknown");
    }

    #[test]
    fn finds_the_mount_for_a_path() {
        let mi = "22 1 8:2 / / rw - ext4 /dev/sda2 rw\n\
                  30 22 8:17 / /data rw - ext4 /dev/sdb1 rw\n\
                  31 22 0:45 / /mnt/ris rw - cifs //storage1/share rw\n\
                  32 22 8:33 / /mnt/my\\040disk rw - xfs /dev/sdc1 rw\n";
        assert_eq!(linux_mount_for("/data/run1.arrow", mi), Some(("ext4".into(), "/dev/sdb1".into())));
        assert_eq!(linux_mount_for("/database/x", mi), Some(("ext4".into(), "/dev/sda2".into())));
        assert_eq!(linux_mount_for("/mnt/ris/a.arrow", mi), Some(("cifs".into(), "//storage1/share".into())));
        assert_eq!(linux_mount_for("/mnt/my disk/a", mi), Some(("xfs".into(), "/dev/sdc1".into())));
    }

    #[test]
    fn missing_paths_take_their_nearest_existing_folder_and_volumes_are_probed_once() {
        let dir = std::env::temp_dir();
        let missing = dir.join("pioneer_not_created_yet").join("results");
        let r = kinds(&[missing.to_string_lossy().into_owned(), dir.to_string_lossy().into_owned()]);
        assert_eq!(r[0].kind, r[1].kind);
        assert_eq!(kinds(&[String::new()])[0].kind, "unknown");
        let dir = std::env::temp_dir();
        let r = kinds(&[dir.to_string_lossy().into_owned(), dir.to_string_lossy().into_owned()]);
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].kind, r[1].kind);
    }
}
