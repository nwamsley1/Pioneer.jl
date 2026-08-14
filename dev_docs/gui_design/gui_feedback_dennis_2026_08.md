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

## Group 2 — MBR parameter — [x] done

- [x] **#4** Surface `global.match_between_runs` as a toggle.

Labelled `MBR` with `Match between runs` beneath, defaulting to on, in
Confidence & output for now — it moves to Common with the rest in Group 3.

Touched: `types.ts` (field + default), `config.ts` (`SEARCH_OWNED_PATHS`,
`buildSearchJsonBase`, `searchConfigToState`), one `ToggleRow` in
`SearchDiaForm.tsx`. `onToggle` in `App.tsx` is generic, so no handler work.

**Default is `true`.** `IntegrateChromatogramsSearch.jl:169-171` falls back to
`true` when the key is absent, and the GUI never emitted it before — so every
GUI run to date has had MBR on. Defaulting the toggle to `false` would have
silently changed behaviour for existing users.

The key is now emitted on every run even though the value matches Pioneer's
fallback: the config is the record of what a run was told to do, and MBR moves
the result enough that "it was left unset" is not a useful answer later. On
load, an absent key leaves the toggle as it stands rather than reading as off.

## Group 3 — Section reorganization — [x] done

- [x] **#5** Rename and regroup: COMMON PARAMETERS / ADVANCED PARAMETERS / OUTPUT.
- [x] **#6** Contents:
      - Common: q-value threshold, run-to-run normalization, MBR
      - Advanced: fragment isotopes, initial NCE, min peptides
      - Output: write CSV, write decoys, delete temp, debug logging

Final layout, after two rounds of review at the running app:

    Files and name        (full width)
    Common parameters  |  Output        (equal heights, side by side)
    Advanced parameters   (full width, beneath both)

The first tried arrangement — tallest card left, other two stacked right — left
about 150px of dead space under the left card, since one card cannot balance
two. Putting Output beside Common and dropping Advanced full-width fixes it:
the two in the row differ by only a row or so, and Advanced's three fields get
the full width instead of ~87px per column, so their labels stop wrapping.
`alignItems: stretch` makes the pair equal height; when the row wraps at narrow
widths each card is alone on its line and returns to its natural height.

Also folded in:

- Dropped "Defaults suit most experiments" from Advanced.
- Dropped "Left blank, a name is generated for you." from Job name. The whole
  line is now conditional, so there is no empty paragraph holding space — the
  "Optional" placeholder already carries the meaning.
- Moved Job name into the first card. It is shared with Convert, Download and
  BuildSpecLib, which have no equivalent card, so it became `JobNameField`:
  inline for SearchDIA, in a card of its own for the rest.
- Renamed "Essentials" to "Files and name" — the old name said nothing about
  what was in the card.
- Dropped "· all settings on one page" from the SearchDIA and BuildSpecLib
  page subtitles. It described the layout rather than the command, and the
  other two subtitles carry no such suffix.

**Still open:** the Results folder lives in "Files and name" while a card called
"Output" holds write-CSV / decoys / temp, so someone looking for where results
are written may check "Output" first. Suggested fix is renaming that card to
"Output options"; not applied.

## Group 4 — Real modification names — [x] done

- [x] **#1** Show UNIMOD names rather than `Unimod:NN` in the library summary.

`unimodDisplay()` in `fasta.ts` rewrites one `library_info` entry, swapping the
accession for the name and leaving the site suffix alone. Used by
`LibrarySummary` for both fixed and variable mods, and by the Prosit
variable-mods warning banner in `SearchDiaForm`, which printed the same raw
accessions in prose.

Done at the display layer rather than by restructuring the Rust return type:
the string is a display artifact either way, and changing it would have touched
the command's serde struct and the four `paths.rs` tests asserting on those
exact strings, for nothing a user sees.

Coverage is the existing `KOINA_MODS` catalogue, which `unimodLabel` already
scans across every model. That covers Altimeter's only two modifications —
Carbamidomethyl (4) and Oxidation (35) — plus the Prosit PTM set. An accession
outside it is left exactly as it came: opaque beats wrong, and the number is
still something the reader can look up.

Verified against the real formats:

    "Unimod:4 on C"      -> "Carbamidomethyl on C"
    "Unimod:35 on M"     -> "Oxidation on M"
    "Unimod:21 on [ST]"  -> "Phospho on [ST]"     (regex site preserved)
    "Unimod:4"           -> "Carbamidomethyl"     (no site recorded)
    "Unimod:9999 on X"   -> unchanged             (unknown accession)
    "Acetyl on K"        -> unchanged             (already a name)

## Digest preview responds to the digestion settings — [x] done

Not from Dennis's list; raised by Nathan while working through it.

The preview under the enzyme picker showed raw cut fragments and ignored
specificity, length limits and missed cleavages, so changing any of them left
it unmoved. `previewDigest` is now a faithful port of `digest_sequence`
(`fasta_digest.jl`) and lists the peptides a library would actually carry.

It leads with the count, because that is the part that moves. On the sample
sequence at the default 7–40 residues and one missed cleavage:

    full    3 peptides
    semi    54
    semi-n  22
    semi-c  35

The order-of-magnitude jump on switching to semi is the most consequential
thing about that choice and was otherwise invisible until a build took all
night. Only the first four peptides are listed, then "+N more".

`cleavageNote` used to judge the rule by counting what `previewDigest`
returned. That no longer works — a sound rule can yield zero peptides under
tight length limits, which is a different complaint — so cut sites are now
exposed separately as `cleavageSites` and the validation uses those.

Verified against Julia: 12 cases (4 specificities x 3 length/missed-cleavage
settings) match peptide-for-peptide, via `Pioneer.digest_sequence` on the same
sample sequence.

## Reject fixed/variable modifications on the same residue — [x] done

Not from Dennis's list; found by Nathan while clicking around.

The form let you set Oxidation as fixed on M while it was already variable on M.
Pioneer rejects that config outright (`check_mod_site_conflicts`), so the build
failed only after being launched. `modSiteConflict` in `validate.ts` ports that
rule and blocks the run inline, scrolling to the modification table the way the
other mod validations do.

The rule is about the *residue*, not the modification: a fixed carbamidomethyl
on C conflicts with a variable oxidation on C, because a fixed modification
occupies every matching residue and the variable one would stack on top of it.
`modPatternResidues` ports `mod_pattern_residues`, so a class like `[ST]`
resolves and a multi-residue rule like `K[^P]` matches no single residue and is
left to another check.

Checked, including the cases that must stay allowed:

    fixed Ox/M        + variable Ox/M          BLOCK  (the reported case)
    fixed Carbam/C    + variable Ox/M          allow  (the default setup)
    fixed Carbam/C    + variable Ox/C          BLOCK  (different mod, same site)
    fixed Phospho/[ST]+ variable Ox/[TY]       BLOCK  (reports only T)
    fixed K[^P]       + variable Ox/K          allow  (matches no single residue)

Both rows involved turn red as the clash is created, not on Run:
`conflictingResidues` is recomputed each keystroke and both tables receive the
same set, since either row can resolve it and marking one would point at the
wrong half as often as not. Red is reserved for this; the existing amber still
means "this model does not accept that modification". The run-block stays as
the backstop.

Note the guard the form already had was one-directional: the variable table
gets `occupied` from the fixed mods, so a variable mod cannot take a fixed
residue, but the fixed table gets an empty set — which is how a fixed mod could
be added onto a variable residue.

## Drag to reorder the queue — [x] done

Not from Dennis's list; raised by Nathan.

The scheduler takes the first `queued` job in the `jobs` array and the sidebar
renders the queue in that same order, so reordering the array reorders
execution — no separate ordering to store. Queued jobs do not survive a restart
(they are re-read as `interrupted`, deliberately), so this is session state
only and needs nothing from `history.rs`.

`moveQueuedJob` in `lib/queue.ts` holds the running and finished jobs' slots
fixed and rewrites only the queued entries into the remaining ones, so a
running job keeps the row it is displayed at. Six cases checked, including a
running job sitting between two queued ones and a drop aimed at a running row.

Affordance: a six-dot grip handle at the row's leading edge, the convention in
Notion, Linear, Jira and Trello. Held at 45% opacity rather than revealed on
hover — a handle nobody can see is a feature nobody finds. `cursor: grab` /
`grabbing`, the carried row fades to 40%, and the row under the pointer takes
an accent line along the edge the drop will land on.

Handles appear only when the sidebar is expanded and more than one job is
queued. A dedicated handle rather than a draggable row, because the row already
opens the job on click and a handle separates the two without a press-and-hold
heuristic.

**Not done:** HTML5 drag-and-drop is not keyboard accessible. `alt+up/down` on a
focused queue row would give parity, roughly ten lines.

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
