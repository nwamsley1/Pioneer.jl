# Pioneer website

The public site for Pioneer: home page, downloads, license and citation. It is
a static [Astro](https://astro.build) site styled with Tailwind CSS, served at
https://nwamsley1.github.io/Pioneer.jl/ from this repository's `gh-pages`
branch.

The documentation and regression reports are **not** built here. Documenter
publishes the docs to `docs/` on `gh-pages`, and the regression workflows
publish the reports to `reports/`. This site links to them; their URLs live in
`src/site.config.ts`.

## Pages

| Page | Source | Data |
|---|---|---|
| `/` | `src/pages/index.astro` | Latest release version from GitHub |
| `/download/` | `src/pages/download.astro` | GitHub Releases API: installers, sizes, release notes, earlier versions |
| `/license/` | `src/pages/license.astro` | `LICENSE` from the latest release tag, falling back to `src/data/LICENSE.fallback.txt` |
| `/cite/` | `src/pages/cite.astro` | `paper` in `src/site.config.ts` |

Installers are matched to platforms by file name: `.pkg` with `arm64` or
`x64`, `.exe`/`.msi`, and `.deb`. If the installer naming in the build
workflows changes, update `classify()` in `src/lib/releases.ts`.

## Local development

```bash
npm ci
npm run dev        # http://localhost:4321/Pioneer.jl/
npm run build      # static output in dist/
npx astro preview  # serve dist/ at http://localhost:4321/Pioneer.jl/
```

Every page is served under `/Pioneer.jl/` (`base` in `astro.config.mjs`), and
internal links use `import.meta.env.BASE_URL`, so they work both locally and on
GitHub Pages.

Set `GITHUB_TOKEN` to avoid the API's anonymous rate limit. Without network
access the build uses `src/data/releases.sample.json` and the Download page
shows a banner saying the data is a sample. To build from a saved copy of the
API response instead:

```bash
curl -sS "https://api.github.com/repos/nwamsley1/Pioneer.jl/releases?per_page=30" -o .cache/releases.json
RELEASES_FILE=.cache/releases.json npm run build
```

## Deployment

`.github/workflows/website.yml` builds the site and commits `dist/` to the root
of the `gh-pages` branch. It runs:

- on every push to `develop` that changes `website/`;
- when a release is published (`release.yml` starts it after "Publish Release");
- daily, as a fallback for picking up new releases;
- on demand from the Actions tab.

The build sets `REQUIRE_RELEASES=true`, so a failed GitHub API request fails the
deploy rather than publishing sample data.

The workflow records the files it wrote in `.website-manifest` at the root of
`gh-pages`. On the next deploy it removes only those files before copying the
new build, and it never touches `docs/`, `dev/` or `reports/`.

Static assets go to `site-assets/` rather than Astro's default `_astro/`:
`gh-pages` is processed by Jekyll, which hides folders whose names start with
an underscore.

## Custom domain

To move the site to its own domain later:

1. In the repository's **Settings → Pages → Custom domain**, enter the domain
   and enable **Enforce HTTPS** once the certificate is issued, then point the
   domain's DNS at GitHub Pages.
2. Set `base: "/"` and the new `site` in `astro.config.mjs`.
3. Set `url` in `src/site.config.ts` to the new domain.

The custom domain applies to the whole `gh-pages` site, so the docs and reports
move with it; update `docsUrl`, `installGuideUrl` and `reportsUrl` in
`src/site.config.ts` to match.

## Before launch

- Confirm with WashU OTM that the wording on the home and License pages is
  acceptable, and that the site may describe the software before the startup
  license is signed.
- Don't add Washington University names or logos beyond the plain-text
  attribution without permission from WashU.
- Check the citation in `src/site.config.ts` (author list and DOI) against the
  current preprint.
