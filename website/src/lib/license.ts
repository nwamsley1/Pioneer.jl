import { site } from "../site.config";
import fallback from "../data/LICENSE.fallback.txt?raw";

// The License page shows the LICENSE file from the tag of the latest release,
// so it always matches what people download. If that can't be fetched, the
// copy bundled with the site is used.
export async function getLicenseText(tag: string | undefined): Promise<string> {
  if (tag) {
    try {
      const res = await fetch(`https://raw.githubusercontent.com/${site.repo}/${tag}/LICENSE`);
      if (res.ok) return (await res.text()).trim();
    } catch {
      // fall through to the bundled copy
    }
  }
  return fallback.trim();
}

// Splits the license into paragraphs, dropping the blank spacer lines.
export function paragraphs(text: string): string[] {
  return text
    .split(/\n\s*\n/)
    .map((p) => p.replace(/\s+/g, " ").trim())
    .filter((p) => p.length > 0);
}
