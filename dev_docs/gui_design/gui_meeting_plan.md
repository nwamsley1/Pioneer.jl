---
title: "Pioneer Console — plan from the GUI meeting"
subtitle: "Twenty-one notes, grouped, costed and sequenced · 2026-08-12"
geometry: margin=1in
fontsize: 11pt
---

# How to read this

Every item below comes from the meeting notes. Each one says what it means,
**what the code does today** (checked, not assumed), and what the work is.

Effort is **S** (an hour or two), **M** (half a day to a day), **L** (multi-day,
or needs a decision first).

Three notes are already partly built, two are bugs with a known cause, one is
not a GUI change at all, and two need a decision from you before they can be
costed. Those are called out rather than folded in.

\newpage

# A. Chrome and layout

## A1. Collapse icon smaller, moved to the top — S

Today the collapse control is a full-width pill at the *bottom* of the sidebar
with a chevron and the word "Collapse", `padding: 12`, `font 600 13.5px`. It is
the heaviest element in the sidebar for something pressed once a session.

Move it to the top bar of the sidebar as an icon-only button, drop the label
and the fill. This pairs with A2 — both free space at the top of the sidebar.

## A2. Remove the large logo — S

`PioneerLogo` sits in its own block with a bottom border, and since the
title-bar change it also carries the traffic-light clearance. Removing it
reclaims roughly 90px at the top of the sidebar.

Worth deciding what, if anything, replaces it. Options: nothing (the window
title and the Dock icon already identify the app), a small wordmark on the same
row as the collapse button, or keep the peak mark at a much smaller size. A1
and A2 should be designed together as one header row.

\newpage

# B. Naming and wording

These are all one-line changes to labels and hints. Grouped because they should
be reviewed together for tone, and shipped as one commit.

## B1. "Calibration file" — S

The name says *how* it is used rather than *what to give it*. Its actual job is
to supply a reference run from which the fragment and precursor m/z ranges are
detected (`auto_detect_frag_bounds`).

Suggested: **Reference MS file**, hint *"Any run from this experiment — used to
detect the fragment and precursor m/z range."*

## B2. "NCE" — S

The field is a starting guess that Pioneer refines, not a fixed instrument
setting. Suggested: **NCE (starting estimate)**, hint noting it is refined
during the search.

## B3. "Isotopes" — S

Ambiguous between precursor and fragment isotopes. The parameter is fragment
isotopes; the label should say so.

## B4. Debug logging level — S/M

A select for the log verbosity Pioneer emits. **Open question:** what levels
does Pioneer actually accept, and is it a CLI flag, an env var, or a config
key? The GUI work is trivial once that is settled; the answer decides whether
this is S or M.

\newpage

# C. Paths and file selection

The largest cluster, and the one with the most existing machinery to build on.

## C1. Library output accepts a folder that does not exist — S

Today `browseLibPath` calls `pickFolder`, so the native dialog can only return a
directory that already exists; the `.poin` name is appended afterwards. That
means you cannot name a new library in the place you want it.

Switch to a save dialog seeded with a default filename, and keep appending
`.poin` when it is missing.

## C2. Configurable default browse location — S (partly built)

**Already done:** `backend.ts` keeps `pioneerConsole.lastDir` in localStorage
and passes it as the dialog's `defaultPath`, so pickers reopen where you last
were, across sessions as well as within one.

**Not done:** there is no configured default and no home-directory fallback —
with no stored value the dialog opens wherever the OS decides. Add a setting
for a preferred root, and fall back to it, then to home, before letting the OS
choose.

## C3. Results folder derived from the MS data folder — S

A button beside *Results folder* that fills it from the MS data folder plus a
suffix. **Open question:** what suffix — a fixed `_results`, the date, or the
run name the GUI already generates? The run name would make results
self-describing but changes on every run, so it cannot be filled in ahead of
time without also pinning the name.

## C4. Remember the last spectral library — S

The library path is the least-changing field in a SearchDIA run and the most
annoying to re-pick. The draft is already persisted under
`pioneerConsole.v2`, so the gap is that a fresh install starts empty rather
than that the value is lost. Cheapest version: keep a most-recently-used list
of libraries and offer it as a dropdown on the field.

## C5. Drag a folder onto the browse bar — M

Tauri v2 gives file-drop events per window; the work is a drop target on each
path field plus the hover affordance, and rejecting drops of the wrong kind
(a file where a folder is wanted). Do it once as a shared component so all
three forms get it.

## C6. Pick a subset of files in the MS data folder — L

Today the field is a folder and Pioneer takes everything in it. Selecting a
subset needs a file list in the UI with checkboxes, and — the harder half — a
way to express the subset to Pioneer. **Open question:** does SearchDIA accept
an explicit file list, or only a directory? If only a directory, the options
are a temporary folder of symlinks or a change on the Julia side. That answer
decides whether this is M or L.

\newpage

# D. Library awareness

## D1. Show metadata from the library — M

Once a `.poin` is selected, read and display what it is: the prediction model,
digest parameters, modifications, precursor count. This turns the path field
from a string into something you can sanity-check before a two-hour run.

The library is a directory with an Arrow table plus config, so this is a Rust
side read of the library's own `config.json` (the inspector already detects
`.poin` layout via `PION_MARKERS`) and a small summary panel under the field.

Note the wording in the meeting was "metadata from the library *name*" — if
what you meant is parsing the filename rather than reading the library, that is
much cheaper but much less reliable. Reading the config is the version worth
building.

## D2. DownloadSpecLib as a fourth workflow tab — L

A new command alongside ConvertRAW / BuildSpecLib / SearchDIA for fetching a
prebuilt library.

This is the largest item in the notes and the only one that needs backend that
does not exist yet. **Open questions:** what is the source (a Pioneer-hosted
index, Zenodo, something else), is there an index to list from or only known
URLs, and does the download need resume and checksum verification. Until those
are settled it cannot be costed beyond "L".

The GUI half is well understood — it is the same form-and-run shape as the
other three tabs — so the risk is entirely in the source and transport.

\newpage

# E. Run lifecycle

## E1. A console window opens on Windows every run — S, real bug

`runner.rs` spawns with a plain `Command::new`. On Windows, a GUI process
starting a console subsystem executable gets a console window unless
`CREATE_NO_WINDOW` (`0x08000000`) is set through
`std::os::windows::process::CommandExt::creation_flags`.

The `taskkill` call in the cancel path needs the same treatment, or cancelling
flashes a second window.

Cheap and self-contained, but only verifiable on a Windows build.

## E2. Open the results folder when a run completes — S

Needs `tauri-plugin-opener`; only `tauri-plugin-dialog` is in `Cargo.toml`
today. Add the plugin, then a button on the completed-run view and in the log
drawer's header.

Worth pairing with the existing check for whether the folder still exists — the
drawer already reports when a restored run's output has been moved or deleted,
and the button should be hidden or disabled in that case rather than opening
nothing.

## E3. Fewer default threads for BuildSpecLib — S, but answer the question first

**Open question, and it should be measured rather than guessed:** how much does
thread count actually affect a library build? The phases are Koina requests
(network-bound, governed by `max_koina_requests`, currently 24) and local
fragment prediction and indexing (CPU-bound).

The 7.2-minute Prosit build I ran earlier spent 300s of its 431s in Fragment
Prediction, which suggests threads do matter there. One build at 4 threads
versus one at 15 would settle it in under half an hour, and the answer is worth
knowing regardless of what the default becomes.

\newpage

# F. State and persistence

## F1. Queue and history numbering — S, needs one clarification

The badge renders `{idx + 1}` over each *filtered* list, so queue and history
number independently and both start at 1. Deleting an entry renumbers
everything below it.

**Which behaviour did you want?** The note reads two ways: either the numbers
should be stable per run (an identity that survives deletion), or they are
currently failing to renumber. From the code, positional renumbering is what it
does. If you want stable identity, that is a per-run counter stored with the
run, which also makes F2 and F3 easier.

## F2. Move history and queue to SQLite — M/L

Today: `localStorage`, with history capped at 100 runs and parameters only —
no log output, because a single SearchDIA log is megabytes against a
low-single-digit-megabyte budget shared with the form draft.

SQLite removes the cap, makes F3 a query instead of a scan, and would let logs
be kept. Needs `tauri-plugin-sql` or `rusqlite`, a schema, and a migration that
imports whatever is in `localStorage` so existing history is not dropped.

Do this before F3 and F4 rather than after — both are much cheaper on top of a
real store, and doing them first means writing them twice.

## F3. Search the run history — S once F2 lands, M before it

A filter over run name, command, target path and date.

## F4. Persist the queue across restarts — M

Two behaviours are possible and they are not the same: re-queue the pending
runs on launch, or move them to history as "interrupted". Re-queuing risks
starting a long job the moment the app opens; the note suggests you are open to
the simpler option. My recommendation is to record them as interrupted and
offer a one-click re-run, which is safe by default and loses nothing.

\newpage

# G. Not a GUI change

## G1. LightGBM slower on Windows than macOS — investigation

This belongs to Pioneer, not the console, but it should be tracked because it
shapes what users experience.

Worth establishing first whether it is real and by how much — same library,
same data, wall-clock on both. If it holds up, the usual suspects are the
threading backend (OpenMP on Windows versus the macOS build), the compiler
flags the Windows binary ships with, and whether the packaged sysimage's
`cpu_target` reaches the same instruction set on both. The portable-packaging
work changed exactly that, so a measurement taken before it will not answer the
question.

\newpage

# Suggested sequencing

**First — the cheap wins and the real bugs.** All small, all independent, all
visible immediately:

- E1 Windows console window (a bug, and it makes the app look broken)
- C1 library output save dialog (a bug, blocks a normal workflow)
- B1–B3 the renames
- A1 + A2 sidebar header, designed as one change
- E2 open results folder

**Second — the path work**, which shares a component and should not be spread
out: C2 configured default, C3 derived results folder, C4 remembered library,
C5 drag and drop.

**Third — F2 SQLite**, then F1, F3 and F4 on top of it. Doing F3 or F4 before
the store means writing them twice.

**Fourth — the two large ones**, once their open questions have answers: D1
library metadata, C6 file subset, D2 DownloadSpecLib.

**In parallel, whenever** — E3 and G1 are measurements, not features. Both can
be run by someone else while the above proceeds, and both change decisions
rather than code.

# Open questions, collected

1. **B4** — what logging levels does Pioneer accept, and how are they set?
2. **C3** — what suffix should the derived results folder use?
3. **C6** — can SearchDIA take an explicit file list, or only a directory?
4. **D1** — read the library's config, or parse its filename?
5. **D2** — what is the download source, and is there an index?
6. **F1** — stable per-run numbers, or positional?
7. **F4** — re-queue interrupted runs, or record them as interrupted?
