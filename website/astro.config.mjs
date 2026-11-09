import { defineConfig } from "astro/config";
import tailwindcss from "@tailwindcss/vite";

export default defineConfig({
  site: "https://pioneeromics.bio",
  // Served at the root of the custom domain set on the Pioneer.jl repository's
  // GitHub Pages site.
  base: "/",
  trailingSlash: "ignore",
  build: {
    // The gh-pages branch is processed by Jekyll, which hides folders that
    // start with an underscore, so Astro's default "_astro" would 404.
    assets: "site-assets",
  },
  vite: {
    plugins: [tailwindcss()],
  },
});
