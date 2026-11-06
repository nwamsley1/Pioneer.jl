// Everything that is likely to change lives here, so moving the docs to a
// subdomain or renaming the repository is a one-line edit.

export const site = {
  name: "Pioneer",
  // Production URL of this site. Used for canonical links and the sitemap.
  url: "https://nwamsley1.github.io",

  // GitHub repository whose releases appear on the Download page.
  repo: "nwamsley1/Pioneer.jl",
  // Releases older than this are hidden everywhere on the site.
  minVersion: "0.12.0",

  // Existing Documenter build and regression reports. Point these at
  // docs.pioneer.bio / reports.pioneer.bio once those subdomains exist.
  docsUrl: "https://nwamsley1.github.io/Pioneer.jl/docs/stable/",
  installGuideUrl:
    "https://nwamsley1.github.io/Pioneer.jl/docs/stable/user_guide/installation/",
  reportsUrl: "https://nwamsley1.github.io/Pioneer.jl/reports/",

  dockerImage: "dennisgoldfarb/pioneer",

  // Commercial licensing goes through WashU OTM until the startup license
  // is in place.
  commercialContact: "otm@wustl.edu",
  contactEmails: ["n.t.wamsley@wustl.edu", "d.goldfarb@wustl.edu"],

  paper: {
    title:
      "Pioneer and Altimeter: Fast Analysis of DIA Proteomics Data Optimized for Narrow Isolation Windows",
    authors: ["Wamsley, N. T.", "Wilkerson, E. M.", "Major, M.", "Goldfarb, D."],
    bibtexAuthors: "Wamsley, N. T. and Wilkerson, E. M. and Major, M. and Goldfarb, D.",
    venue: "bioRxiv",
    year: 2026,
    doi: "10.64898/2026.02.16.706201",
    url: "https://www.biorxiv.org/content/10.64898/2026.02.16.706201v3",
  },
} as const;
