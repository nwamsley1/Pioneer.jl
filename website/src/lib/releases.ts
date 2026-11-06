import { marked } from "marked";
import { site } from "../site.config";
import sample from "../data/releases.sample.json";

// Release data is read from the GitHub Releases API once per build.
//
// In CI (REQUIRE_RELEASES=true) a failed request fails the build, so a broken
// Download page is never deployed. Locally, if the API can't be reached, the
// bundled sample is used and the page shows a banner saying so.

export type Platform = "macos-arm64" | "macos-x64" | "windows" | "linux";

export interface Asset {
  name: string;
  url: string;
  size: number;
  platform: Platform;
}

export interface Release {
  tag: string;
  version: string;
  name: string;
  date: Date;
  htmlUrl: string;
  // GitHub's automatic source archives for the release tag.
  sourceZip: string;
  sourceTarGz: string;
  sourceTree: string;
  notesHtml: string;
  assets: Asset[];
}

export interface ReleaseData {
  latest: Release | null;
  previous: Release[];
  isSample: boolean;
}

interface GhAsset {
  name: string;
  browser_download_url: string;
  size: number;
}

interface GhRelease {
  tag_name: string;
  name: string | null;
  published_at: string | null;
  html_url: string;
  body: string | null;
  draft: boolean;
  prerelease: boolean;
  assets: GhAsset[];
}

export const platformLabels: Record<Platform, string> = {
  "macos-arm64": "macOS, Apple silicon",
  "macos-x64": "macOS, Intel",
  windows: "Windows",
  linux: "Linux (Debian, Ubuntu)",
};

export const platformOrder: Platform[] = ["macos-arm64", "macos-x64", "windows", "linux"];

function classify(name: string): Platform | null {
  const n = name.toLowerCase();
  if (n.endsWith(".pkg") || n.endsWith(".dmg")) {
    if (/(arm64|aarch64|apple)/.test(n)) return "macos-arm64";
    if (/(x64|x86_64|intel)/.test(n)) return "macos-x64";
    return "macos-arm64";
  }
  if (n.endsWith(".exe") || n.endsWith(".msi")) return "windows";
  if (n.endsWith(".deb")) return "linux";
  return null;
}

function toRelease(r: GhRelease): Release {
  const assets: Asset[] = [];
  for (const a of r.assets) {
    const platform = classify(a.name);
    if (platform) {
      assets.push({ name: a.name, url: a.browser_download_url, size: a.size, platform });
    }
  }
  // Within a platform, list the bundle installer (.exe) before a bare .msi.
  assets.sort((x, y) => (x.name.endsWith(".exe") ? -1 : 0) - (y.name.endsWith(".exe") ? -1 : 0));
  return {
    tag: r.tag_name,
    version: r.tag_name.replace(/^v/, ""),
    name: r.name || r.tag_name,
    date: new Date(r.published_at ?? Date.now()),
    htmlUrl: r.html_url,
    sourceZip: `https://github.com/${site.repo}/archive/refs/tags/${r.tag_name}.zip`,
    sourceTarGz: `https://github.com/${site.repo}/archive/refs/tags/${r.tag_name}.tar.gz`,
    sourceTree: `https://github.com/${site.repo}/tree/${r.tag_name}`,
    notesHtml: marked.parse(r.body ?? "", { async: false }) as string,
    assets,
  };
}

async function fetchFromGitHub(): Promise<GhRelease[]> {
  // RELEASES_FILE=path/to/releases.json builds from a saved copy of the API
  // response (e.g. `curl .../releases?per_page=30 -o .cache/releases.json`),
  // for builds made where api.github.com can't be reached.
  const file = process.env.RELEASES_FILE;
  if (file) {
    const { readFile } = await import("node:fs/promises");
    return JSON.parse(await readFile(file, "utf8")) as GhRelease[];
  }

  const headers: Record<string, string> = {
    Accept: "application/vnd.github+json",
    "X-GitHub-Api-Version": "2022-11-28",
    "User-Agent": "pioneer-site-build",
  };
  const token = process.env.GITHUB_TOKEN;
  if (token) headers.Authorization = `Bearer ${token}`;

  const res = await fetch(`https://api.github.com/repos/${site.repo}/releases?per_page=30`, {
    headers,
  });
  if (!res.ok) throw new Error(`GitHub API returned ${res.status} for ${site.repo} releases`);
  const data = await res.json();
  if (!Array.isArray(data)) throw new Error("GitHub API returned an unexpected response");
  return data as GhRelease[];
}

// Only official numbered releases (v1.2.3 or 1.2.3) appear on the site. Dev
// builds, CI fixture releases and anything else on the Releases page are skipped.
const OFFICIAL_TAG = /^v?\d+\.\d+\.\d+$/;

// Compares dotted version numbers ("0.12.0" vs "0.10.1").
function compareVersions(a: string, b: string): number {
  const pa = a.split(/[.-]/).map((n) => parseInt(n, 10) || 0);
  const pb = b.split(/[.-]/).map((n) => parseInt(n, 10) || 0);
  for (let i = 0; i < Math.max(pa.length, pb.length); i++) {
    const d = (pa[i] ?? 0) - (pb[i] ?? 0);
    if (d !== 0) return d;
  }
  return 0;
}

let cached: Promise<ReleaseData> | undefined;

export function getReleases(): Promise<ReleaseData> {
  cached ??= (async () => {
    let raw: GhRelease[];
    let isSample = false;
    try {
      raw = await fetchFromGitHub();
    } catch (err) {
      if (process.env.REQUIRE_RELEASES === "true") throw err;
      console.warn(`[releases] ${(err as Error).message}; using sample data`);
      raw = sample as GhRelease[];
      isSample = true;
    }
    const published = raw
      .filter((r) => !r.draft && !r.prerelease && OFFICIAL_TAG.test(r.tag_name))
      .map(toRelease)
      .filter((r) => compareVersions(r.version, site.minVersion) >= 0)
      .sort((a, b) => compareVersions(b.version, a.version));
    return { latest: published[0] ?? null, previous: published.slice(1), isSample };
  })();
  return cached;
}

export function formatSize(bytes: number): string {
  if (bytes >= 1e9) return `${(bytes / 1e9).toFixed(1)} GB`;
  return `${Math.round(bytes / 1e6)} MB`;
}

export function formatDate(d: Date): string {
  return d.toLocaleDateString("en-US", { year: "numeric", month: "long", day: "numeric" });
}
