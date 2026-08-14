# GUI feedback (Dennis, 2026-08-14)

Eleven items from Dennis, grouped so each group lands as one commit. Tick a box
when its commit is in. Groups are ordered low-risk first; the two that need a
design decision are last.

Branch: `feat/gui-dennis-feedback`, off `develop` at `6315eec09`.

---

## Group 1 — Copy edits — [x] done

No schema change, no new state.

- [x] **#2** Debug logging hint: drop `" — the log file keeps its own detail either way"`.
      `SearchDiaForm.tsx`. The clause was describing that `debug_file_level`
      stays at 1 regardless of the toggle — true, but nothing the user can act on.
- [x] **#3** Run-to-run normalization hint: `"Retention-time dependent cross-run
      intensity scaling"` (adds "intensity"). `SearchDiaForm.tsx`.
- [x] **#7** Min peptides info icon: add `info` to `NUM_SPECS.minPeptides` in
      `validate.ts`; `NumField` already renders `<InfoDot>` whenever a spec has one.
      Semantics confirmed: `filter_by_min_peptides` is applied inside the per-file
      loop (`protein_inference_pipeline.jl:1929`), so the count is peptides **in
      that run**, not across the experiment.

## Group 2 — MBR parameter — [ ] todo

- [ ] **#4** Surface `global.match_between_runs` as a toggle.

Touch points: `types.ts` (field + default), `config.ts` (`SEARCH_OWNED_PATHS`,
`buildSearchJsonBase`, `searchConfigToState`), one `ToggleRow` in
`SearchDiaForm.tsx`. `onToggle` in `App.tsx` is generic, so no handler work.

**Default must be `true`.** `IntegrateChromatogramsSearch.jl:169-171` falls back
to `true` when the key is absent, and the GUI never emits it today — so every
GUI run currently has MBR on. Defaulting the toggle to `false` would silently
change behaviour for existing users.

## Group 3 — Section reorganization — [ ] todo

- [ ] **#5** Rename and regroup: COMMON PARAMETERS / ADVANCED PARAMETERS / OUTPUT.
- [ ] **#6** Contents:
      - Common: q-value threshold, run-to-run normalization, MBR
      - Advanced: fragment isotopes, initial NCE, min peptides
      - Output: write CSV, write decoys, delete temp, debug logging

Layout (Nathan): the tallest card takes the left half; the other two stack in the
right column. Build it, then critique from the screenshot.

Constraint to respect: the current two cards use `flex: 1 1 300px`, chosen so
they fit at the 980px minimum window width (see the comment above the row in
`SearchDiaForm.tsx`). A two-column split with a stacked right column keeps that
budget. Advanced's three `NumField`s currently sit in a `1fr 1fr 1fr` grid and
will need checking at half width.

## Group 4 — Real modification names — [ ] todo

- [ ] **#1** Show UNIMOD names rather than `Unimod:NN` in the library summary.

`unimodLabel()` already exists (`fasta.ts`) and maps an id to a name via the
`KOINA_MODS` tables. Rust builds display strings like `"Unimod:35 on M"`
(`paths.rs`, `mod_list`), which `LibrarySummary.tsx` prints verbatim.

Plan: rewrite `Unimod:NN` → name at the display layer, falling back to the
accession when the id is unknown — never guess a name. The same helper applies to
the Prosit variable-mods warning banner in `SearchDiaForm.tsx`.

Open: coverage. `KOINA_MODS` has ~29 entries; an Altimeter library can carry any
accession and there is no UNIMOD table elsewhere in the repo. Either (i) ship the
~29 and leave accessions for the rest, or (ii) generate a fuller id→name table
from `unimod.obo`. Leaning (i), revisit if Dennis hits unnamed mods.

## Group 5 — Manual m/z range — [ ] todo

- [ ] **#9** Let the user set fragment/precursor m/z bounds instead of supplying
      a reference file.

Server-side defaults are in `defaultBuildLibParams.json`:
`auto_detect_frag_bounds: true`, `frag_mz 150–2020`, `prec_mz 390–1010`.
`get_frag_bounds.jl` uses the reference file only when auto-detect is on *and*
the file exists; otherwise it warns and falls back to those numbers.

Work: auto/manual control on the existing "Reference MS file" card, four
`NumField`s behind the manual choice, new `BuildParams` fields, `NUM_SPECS`
entries, and the five keys through `BUILD_OWNED_PATHS` / `buildLibJsonBase` /
the loader. Manual sets `auto_detect_frag_bounds: false` and emits the numbers.

## Group 6 — Rename a run in the queue — [ ] todo

- [ ] **#10** Edit icon beside the run name; clicking makes the text editable.

`Sidebar.tsx` renders `{j.title}` in `renderRow`. Needs an inline input plus a
`renameJob(id, title)` in `App.tsx`, reusing `resolveRunName` for collisions.

**Catch:** `history.rs` upserts with `ON CONFLICT(id) DO UPDATE SET status,
target, snapshot` — `title` is not in the update list, so a rename would vanish
on restart. Add `title = excluded.title` there.
