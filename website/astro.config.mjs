import { defineConfig } from "astro/config";
import tailwindcss from "@tailwindcss/vite";

export default defineConfig({
  site: "https://pioneeromics.bio",
  // Served at the root of the custom domain set on the Pioneer.jl repository's
  // GitHub Pages site.
  base: "/",
  trailingSlash: "ignore",
  // Keep the whitespace between text and inline elements (e.g. "under a <a>license</a>"); compressing it
  // joins words across a line break before a link.
  compressHTML: false,
  build: {
    // The gh-pages branch is processed by Jekyll, which hides folders that
    // start with an underscore, so Astro's default "_astro" would 404.
    assets: "site-assets",
  },
  vite: {
    plugins: [tailwindcss()],
  },
});
