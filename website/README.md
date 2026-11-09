# Pioneer website

The public site for Pioneer: home page, downloads, license and citation. It is
a static [Astro](https://astro.build) site styled with Tailwind CSS, served at
https://pioneeromics.bio/ from this repository's `gh-pages` branch.

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
npm run dev        # http://localhost:4321/
npm run build      # static output in dist/
npx astro preview  # serve dist/ at http://localhost:4321/
```

Internal links use `import.meta.env.BASE_URL` (`base` in `astro.config.mjs`),
so the site can move under a sub-path without editing pages.

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

## Download form

The Download page asks for name, email and institution before showing the
installers, and remembers in the browser that the visitor has filled it in.
The form posts to the Cloudflare Worker in `worker/`, which adds everyone to
the Mailchimp audience with the tag `pioneer-download`: visitors who tick the
release-emails box are subscribed, and the rest are stored as non-subscribed
contacts that receive no campaigns.

Setup, once:

1. In Mailchimp, add an audience text field named "Institution" with the merge
   tag `ORG`, and copy the audience ID into `MAILCHIMP_LIST_ID` in
   `worker/wrangler.toml`.
2. In `worker/`: `npx wrangler login`, `npx wrangler secret put MAILCHIMP_API_KEY`
   (a Mailchimp API key), then `npx wrangler deploy`.
3. Put the Worker's URL in `downloadFormEndpoint` in `src/site.config.ts`.

With `downloadFormEndpoint` empty, the form only logs the details to the
browser console.

## Custom domain

`pioneeromics.bio` is set in the repository's **Settings → Pages → Custom
domain**, and its DNS (on Cloudflare) points at GitHub Pages. The domain
applies to the whole `gh-pages` site, so the docs and reports are served from
it too, and the old `nwamsley1.github.io/Pioneer.jl/` URLs redirect to it.
`public/CNAME` keeps the domain in every website deploy.

## Before launch

- Confirm with WashU OTM that the wording on the home and License pages is
  acceptable, and that the site may describe the software before the startup
  license is signed.
- Don't add Washington University names or logos beyond the plain-text
  attribution without permission from WashU.
- Check the citation in `src/site.config.ts` (author list and DOI) against the
  current preprint.
