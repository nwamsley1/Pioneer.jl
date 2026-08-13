#!/usr/bin/env python3

from __future__ import annotations

import argparse
import datetime as dt
import html
import json
import shutil
import subprocess
import sys
from pathlib import Path

try:
    from zoneinfo import ZoneInfo
except ImportError:  # pragma: no cover - Python < 3.9 fallback
    ZoneInfo = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Publish a regression report to gh-pages and prune stale plot assets."
    )
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--report-path", required=True)
    parser.add_argument("--report-subdir", required=True)
    parser.add_argument("--history-branch", default="gh-pages")
    parser.add_argument("--github-repository")
    parser.add_argument("--github-token")
    parser.add_argument("--repo-url")
    parser.add_argument("--retention-days", type=int, default=30)
    parser.add_argument("--current-timestamp")
    return parser.parse_args()


def run_command(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, check=check, text=True, capture_output=True)


def ensure_success(result: subprocess.CompletedProcess[str], context: str) -> None:
    if result.returncode == 0:
        return
    output = (result.stdout or "") + (result.stderr or "")
    raise RuntimeError(f"{context} failed:\n{output}".rstrip())


def iso_utc_now(current_timestamp: str | None) -> str:
    if current_timestamp:
        return current_timestamp
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_timestamp(value: str | None) -> dt.datetime | None:
    if not value:
        return None
    normalized = value.strip()
    if normalized.endswith("Z"):
        normalized = normalized[:-1] + "+00:00"
    try:
        parsed = dt.datetime.fromisoformat(normalized)
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    return parsed.astimezone(dt.timezone.utc)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def derive_repo_url(args: argparse.Namespace) -> str:
    if args.repo_url:
        return args.repo_url
    if not args.github_repository or not args.github_token:
        raise ValueError("Either --repo-url or both --github-repository and --github-token are required.")
    return f"https://x-access-token:{args.github_token}@github.com/{args.github_repository}.git"


def ensure_checkout(repo_dir: Path, repo_url: str, history_branch: str) -> None:
    if (repo_dir / ".git").is_dir():
        fetch_result = run_command(
            ["git", "-C", str(repo_dir), "fetch", "origin", history_branch],
            check=False,
        )
        if fetch_result.returncode != 0:
            combined = (fetch_result.stdout or "") + (fetch_result.stderr or "")
            print(f"warning: fetch failed for {history_branch}: {combined.strip()}", file=sys.stderr)

        checkout_result = run_command(
            ["git", "-C", str(repo_dir), "checkout", history_branch],
            check=False,
        )
        if checkout_result.returncode != 0:
            ensure_success(
                run_command(["git", "-C", str(repo_dir), "checkout", "-b", history_branch]),
                f"checkout branch {history_branch}",
            )

        remote_exists = run_command(
            ["git", "-C", str(repo_dir), "rev-parse", "--verify", f"origin/{history_branch}"],
            check=False,
        )
        if remote_exists.returncode == 0:
            ensure_success(
                run_command(
                    ["git", "-C", str(repo_dir), "reset", "--hard", f"origin/{history_branch}"]
                ),
                f"reset branch {history_branch}",
            )
        return

    if repo_dir.exists():
        shutil.rmtree(repo_dir)

    clone_result = run_command(
        [
            "git",
            "clone",
            "--branch",
            history_branch,
            "--depth",
            "1",
            "--filter=blob:none",
            repo_url,
            str(repo_dir),
        ],
        check=False,
    )
    if clone_result.returncode == 0:
        return

    repo_dir.mkdir(parents=True, exist_ok=True)
    ensure_success(run_command(["git", "init", str(repo_dir)]), "initialize git repository")
    ensure_success(
        run_command(["git", "-C", str(repo_dir), "checkout", "-b", history_branch]),
        f"create branch {history_branch}",
    )
    ensure_success(
        run_command(["git", "-C", str(repo_dir), "remote", "add", "origin", repo_url]),
        "add origin remote",
    )


def publish_report(
    repo_dir: Path,
    report_path: Path,
    report_subdir: Path,
    report_timestamp: str,
    retention_days: int,
) -> Path:
    report_dir = repo_dir / report_subdir
    report_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(report_path, report_dir / "index.html")

    plots_source = report_path.parent / "fdr_plots"
    plots_dest = report_dir / "fdr_plots"
    if plots_dest.exists():
        shutil.rmtree(plots_dest)
    plots_available = plots_source.is_dir()
    if plots_available:
        shutil.copytree(plots_source, plots_dest)

    write_json(
        report_dir / "meta.json",
        {
            "plots_available": plots_available,
            "plots_pruned": False,
            "plots_retention_days": retention_days,
            "timestamp": report_timestamp,
        },
    )
    return report_dir


def sync_landing_page(run_dir: Path, repo_dir: Path) -> None:
    landing_source = run_dir / "pioneer" / "pages" / "index.html"
    logo_source = run_dir / "pioneer" / "figures" / "PIONEER_LOGO.svg"
    if not landing_source.is_file():
        raise FileNotFoundError(f"Missing landing page at {landing_source}")
    if not logo_source.is_file():
        raise FileNotFoundError(f"Missing logo at {logo_source}")

    assets_dir = repo_dir / "assets"
    assets_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(landing_source, repo_dir / "index.html")
    shutil.copy2(logo_source, assets_dir / "pioneer-logo.svg")


def prune_stale_plots(
    repo_dir: Path,
    current_report_dir: Path,
    report_timestamp: str,
    retention_days: int,
) -> list[Path]:
    reports_root = repo_dir / "reports"
    if not reports_root.is_dir():
        return []

    cutoff = parse_timestamp(report_timestamp)
    if cutoff is None:
        raise ValueError(f"Invalid current timestamp: {report_timestamp}")
    cutoff -= dt.timedelta(days=retention_days)

    pruned: list[Path] = []
    for meta_file in sorted(reports_root.rglob("meta.json")):
        report_dir = meta_file.parent
        if report_dir == current_report_dir:
            continue

        metadata = load_json(meta_file)
        report_time = parse_timestamp(str(metadata.get("timestamp", "")))
        if report_time is None:
            continue

        plots_dir = report_dir / "fdr_plots"
        metadata_changed = False
        if plots_dir.is_dir() and report_time < cutoff:
            shutil.rmtree(plots_dir)
            metadata["plots_available"] = False
            metadata["plots_pruned"] = True
            metadata["plots_pruned_at"] = report_timestamp
            metadata["plots_retention_days"] = retention_days
            write_json(meta_file, metadata)
            pruned.append(report_dir)
            continue

        plots_available = plots_dir.is_dir()
        if metadata.get("plots_available") != plots_available:
            metadata["plots_available"] = plots_available
            metadata_changed = True
        if metadata_changed:
            write_json(meta_file, metadata)

    return pruned


def chicago_timestamp(value: dt.datetime | None) -> str:
    if value is None:
        return ""
    if ZoneInfo is None:
        return value.strftime("%Y-%m-%dT%H:%M:%S")
    try:
        local_value = value.astimezone(ZoneInfo("America/Chicago"))
    except Exception:
        local_value = value
    return local_value.strftime("%Y-%m-%dT%H:%M:%S")


def render_reports_index(repo_dir: Path) -> None:
    reports_root = repo_dir / "reports"
    reports_root.mkdir(parents=True, exist_ok=True)

    entries: list[dict[str, object]] = []
    for meta_file in reports_root.rglob("meta.json"):
        report_dir = meta_file.parent
        relative_dir = report_dir.relative_to(reports_root)
        parts = relative_dir.parts
        if len(parts) < 4:
            continue

        metadata = load_json(meta_file)
        report_time = parse_timestamp(str(metadata.get("timestamp", "")))
        plots_dir = report_dir / "fdr_plots"
        if plots_dir.is_dir():
            plots_status = "Available"
        elif metadata.get("plots_pruned"):
            plots_status = "Archived"
        else:
            plots_status = "Unavailable"

        entries.append(
            {
                "branch": parts[0],
                "sha": parts[1],
                "regression_set": parts[2],
                "run": parts[3],
                "date_cst": chicago_timestamp(report_time),
                "href": relative_dir.as_posix() + "/",
                "plots_status": plots_status,
                "sort_key": report_time or dt.datetime.min.replace(tzinfo=dt.timezone.utc),
            }
        )

    entries.sort(key=lambda entry: entry["sort_key"], reverse=True)

    rows: list[str] = []
    for entry in entries:
        rows.extend(
            [
                "<tr>",
                f"<td>{html.escape(str(entry['date_cst']))}</td>",
                f"<td>{html.escape(str(entry['branch']))}</td>",
                f"<td>{html.escape(str(entry['sha']))}</td>",
                f"<td>{html.escape(str(entry['regression_set']))}</td>",
                f"<td>{html.escape(str(entry['run']))}</td>",
                f"<td>{html.escape(str(entry['plots_status']))}</td>",
                f"<td><a href=\"{html.escape(str(entry['href']))}\">Open report</a></td>",
                "</tr>",
            ]
        )

    html_output = "\n".join(
        [
            "<!DOCTYPE html>",
            "<html lang=\"en\">",
            "<head>",
            "<meta charset=\"UTF-8\">",
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">",
            "<title>Regression Report Index</title>",
            "<style>",
            "body { font-family: Arial, sans-serif; margin: 24px; }",
            "table { border-collapse: collapse; width: 100%; }",
            "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }",
            "th { background-color: #f4f4f4; }",
            "tr:nth-child(even) { background-color: #fafafa; }",
            "</style>",
            "</head>",
            "<body>",
            "<h1>Regression Reports</h1>",
            "<table>",
            "<thead><tr><th>Date (CST)</th><th>Branch</th><th>Commit</th><th>Regression Set</th><th>Run</th><th>Plots</th><th>Report</th></tr></thead>",
            "<tbody>",
            *rows,
            "</tbody>",
            "</table>",
            "</body>",
            "</html>",
        ]
    )
    (reports_root / "index.html").write_text(html_output + "\n", encoding="utf-8")


def commit_and_push(repo_dir: Path, history_branch: str, run_dir: Path) -> bool:
    ensure_success(run_command(["git", "-C", str(repo_dir), "add", "-A"]), "stage gh-pages updates")
    diff_result = run_command(
        ["git", "-C", str(repo_dir), "diff", "--cached", "--quiet"],
        check=False,
    )
    if diff_result.returncode == 0:
        return False
    if diff_result.returncode not in (0, 1):
        ensure_success(diff_result, "check for staged changes")

    ensure_success(
        run_command(
            [
                "git",
                "-C",
                str(repo_dir),
                "-c",
                "user.name=github-actions[bot]",
                "-c",
                "user.email=github-actions[bot]@users.noreply.github.com",
                "commit",
                "-m",
                f"Update Pages history from {run_dir.name}",
            ]
        ),
        "commit gh-pages updates",
    )
    ensure_success(
        run_command(["git", "-C", str(repo_dir), "push", "origin", history_branch]),
        "push gh-pages updates",
    )
    return True


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    report_path = Path(args.report_path).resolve()
    if not report_path.is_file():
        raise FileNotFoundError(f"Missing report HTML at {report_path}")
    if args.retention_days < 0:
        raise ValueError("--retention-days must be non-negative")

    repo_url = derive_repo_url(args)
    repo_dir = run_dir / "gh-pages-site"
    history_branch = args.history_branch
    report_subdir = Path(args.report_subdir)
    report_timestamp = iso_utc_now(args.current_timestamp)

    ensure_checkout(repo_dir, repo_url, history_branch)
    report_dir = publish_report(repo_dir, report_path, report_subdir, report_timestamp, args.retention_days)
    sync_landing_page(run_dir, repo_dir)
    pruned = prune_stale_plots(repo_dir, report_dir, report_timestamp, args.retention_days)
    render_reports_index(repo_dir)
    pushed = commit_and_push(repo_dir, history_branch, run_dir)

    print(f"Published report: {report_subdir.as_posix()}")
    print(f"Pruned plot directories: {len(pruned)}")
    for report_dir in pruned:
        print(f"  - {report_dir.relative_to(repo_dir).as_posix()}")
    print(f"Changes pushed: {'yes' if pushed else 'no'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
