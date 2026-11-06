import { defineConfig } from "astro/config";
import tailwindcss from "@tailwindcss/vite";

export default defineConfig({
  site: "https://nwamsley1.github.io",
  // Served from the Pioneer.jl repository's GitHub Pages site, so every page
  // lives under /Pioneer.jl/. With a custom domain, set base to "/".
  base: "/Pioneer.jl/",
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
