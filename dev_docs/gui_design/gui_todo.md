# Pioneer Console — remaining work

Checklist derived from the meeting plan (`gui_meeting_plan.pdf`) and the design
notes (`gui_design_notes.pdf`). Tick items as they land.

Last updated 2026-08-12 (third pass).

---

## Next up

- [ ] **Clear the 20 seeded runs** when they have served their purpose.
      Deleting them in the UI works. The run counter deliberately will not
      rewind, so real runs continue from 21.

- [ ] **Decide on MS file names in search** *(raised, not decided)*
      Search now covers every path in a run's snapshot, but not the individual
      files inside the MS data folder — those were never recorded. Capturing
      them means `PathInfo` returning names as well as counts, a new column for
      them, and a migration. It only helps runs recorded from then on; existing
      history stays unsearchable by file name.

- [ ] **Decide on precursor count in the library panel** *(raised, not decided)*
      The single most useful number for "is this the big human library or the
      small ecoli one", but it lives in `precursors_table.arrow` and would mean
      an Arrow reader in the Rust side. A real dependency decision rather than
      something to slip in.

- [ ] **ConvertRAW's Advanced disclosure** — SearchDIA's is now always open,
      since collapsed it left a half-width empty box. ConvertRAW is a single
      column where a collapsed section still reads correctly, so it was left
      alone. Worth matching if you want the consistency.

## Small, ready whenever

- [ ] **E3 — threads for BuildSpecLib** *(a measurement, not a feature)*
      Does thread count actually affect a library build? One build at 4 threads
      against one at 15 settles it. The Prosit build measured earlier spent
      300 s of 431 s in Fragment Prediction, so the answer is probably yes —
      which would mean the default should stay high, not drop.

- [ ] **History diff encoding** *(much less pressing since F2)*
      Store only what differs from the defaults. Measured: SearchDIA 528 B →
      278 B (−47%), BuildSpecLib 959 B → 280 B (−71%); the BuildSpecLib saving
      is mostly the four UniProt header regexes, which `presetId: 'uniprot'`
      already implies. `target` also duplicates `results` and can go.
      The cap and the quota cliff are both gone now that history is in SQLite,
      so this is a size optimisation rather than a correctness one.
      **Must** stamp a defaults version alongside, or changing a default later
      silently reinterprets every stored entry.

---

## Parked

- [ ] **D2 — DownloadSpecLib as a fourth tab**
      Waiting on a download location. The GUI half is the same form-and-run
      shape as the other three tabs; all the risk is in the source and
      transport (is there an index to list, does it need resume and checksum
      verification). A mock against a stub source can be built any time it is
      useful for reviewing the tab's shape.

- [ ] **C6 — pick a subset of files in the MS data folder** *(do last)*
      Needs a file list with checkboxes, and — the harder half — a way to
      express the subset to Pioneer. **Open:** does SearchDIA accept an
      explicit file list, or only a directory? If only a directory, the choices
      are a temporary folder of links or a change on the Julia side.
      Last because it is the only item that may require a Pioneer-side change,
      and nothing else should wait behind a cross-repo dependency.

---

## Not GUI work, but blocking release

- [ ] **Push the branch.** `feat/gui-console-prosit-quant` is well ahead of
      origin and everything recent is local only.

- [ ] **Get a build green.** Our GUI and Prosit work has never been through the
      packaged-binary path with Dennis's portability fixes underneath. As of the
      last check his Windows build was still failing on
      `feature/portable-packaging-phase1`; ours inherits that work, so this is
      gated on it.

- [ ] **Verify the two Windows-only fixes on a real build** — neither can be
      checked from macOS:
      - `CREATE_NO_WINDOW`, so a console no longer opens on every run. The
        Windows branch was type-checked against `x86_64-pc-windows-gnu`, but
        whether the window is gone is unverified.
      - The title bar. `titleBarStyle: "Overlay"` is macOS-only and the
        reserved strip is gated on user-agent, so Windows and Linux should show
        ordinary decorations with no dead space. Untested.

- [ ] **`PioneerConverter` is missing from `build/Pioneer_local`**, so the
      ConvertRAW tab cannot be exercised against the local distribution.

---

## Done

Design notes (`gui_design_notes.pdf`) — §1 stepper, §2 label hierarchy,
§3 field width, §4 Browse weight, §5 toggle scale, §6 stacked card, §7 title
bar, §9 swatch affordance. §8 (icon consistency) withdrawn: the premise was
wrong, all three icons are SVGs at matching stroke weight.

Meeting plan — A1 collapse icon, A2 logo removed, B1 Reference MS file,
B2 Initial NCE, B4 debug logging, C1 library save dialog, C2 browse-from-home
and configurable, C4 remember the last library, B3 Fragment isotopes, C5 drag and drop onto
the field, D1 library metadata, E1 Windows console window, E2 open the output
folder, F1 stable history numbers and no cap, F2 SQLite, F3 history search,
F4 interrupted runs.

Also landed, not from either document — the Koina modification picker, the
settings panel behind the gear, the threads field being clearable, the
single-instance guard, Confidence & output beside Advanced, Advanced always
open on SearchDIA, an explanation of what history search matches, and the
removal of the trace mode control.

A third bug found only by running it: the trace mode control's "separated"
value was never one Pioneer accepts — it wants "combined" or "separate", so
choosing it failed the run at IntegrateChromatogramsSearch. Removing the
control removed the bug.

Two bugs found by testing rather than by reading: the window could not be
dragged at all (the overlay title bar needs core:window:allow-start-dragging,
which core:window:default does not grant), and drops did nothing (Tauri types
the drop position as physical but reports logical pixels, so dividing by the
device ratio halved every coordinate).

---

# Review round, 2026-08-13

Raised in Slack after driving the app. Each item records what the code
actually does today, so the work is scoped before it is started rather than
rediscovered. Tick as they land.

## 1. ConvertRAW threading — **done**

**Ask:** no thread picker in the header here; recover the threads-per-file
setter that older versions had in the options menu. Files-at-a-time stays at 1.

- [x] Header thread picker hidden on ConvertRAW.
- [x] `threadsPerFile` restored to the Advanced grid, default 3.
- [x] Args read the form field rather than the sidebar count.

The picker drove `JULIA_NUM_THREADS`, which PioneerConverter never reads — it
is a .NET program — so on this page it was a control that changed nothing.
`concurrentFiles` was deliberately *not* restored: the two knobs multiply, and
one file at a time is the wanted behaviour. Definitions recovered from
`cb014888f^`, the commit that removed them.

## 2. BuildSpecLib thread count — *last, after measurement*

- [ ] Investigate where threading actually helps library building, pick a
      sensible number, then decide what (if anything) to expose.

Deliberately deferred: this is a measurement task, not a UI task.

## 3. Enzyme presets, plus custom regex

- [ ] Presets for the common enzymes, each carrying its cleavage regex.
- [ ] A custom option letting a regex be typed, validated before it can run.

Today the digestion regex is a single free-text field. Semi-enzymatic is out
of scope — Dennis is working on it separately.

Validation needs deciding: a regex that compiles is not necessarily one
Pioneer accepts, so "valid" may mean more than "parses".

## 4. Carbamidomethyl (C) must be a required fixed mod — **done**

- [x] Always present as a fixed mod, on C, for every supported Koina model.
- [x] Not offered as a variable modification at all.
- [x] Its row carries no remove control.

The rule lives in `lib/koinaMods.ts` (`REQUIRED_FIXED_UNIMOD`,
`requiresFixedAlkylation`, `enforceRequiredMods`) rather than in the form, so a
future model that does not want it opts out in one place. The model set is
listed explicitly, so adding a model is a decision rather than an inheritance.

`enforceRequiredMods` runs wherever the lists change from outside the editor —
a model switch, a loaded config — and normalises six cases: absent, present,
wrongly variable, duplicated, fixed on the wrong site, and variable on a PTM
model.

The pin is on the **site**, not the modification: `[CK]` is a valid fixed
pattern, K-carbamidomethyl is freely available as fixed or variable on the PTM
models, and only C is reserved. A fixed row covering K alone gains C rather
than being replaced, so a chosen modification is never silently dropped.

## 5. Queue items need the history items' descriptor — **done**

- [x] Descriptor shown on running rows too.

Both lists share one `renderRow`; the subtext was the else-branch of
`running`, so the progress bar replaced it. A running row was therefore the
one place the command and status were not written down, while being the row
most worth identifying. The bar now sits beneath the descriptor instead.

## 6. File count in the subtext

- [ ] Show the number of MS files under queue and history rows.

`PathInfo` already carries `entry_count`, but a finished run's count is not
stored — `StoredRun` would need the number recorded at enqueue time, since the
folder may have changed by the time the row is drawn.

## 7. Row descriptor format

- [ ] `SearchDIA · N files · Completed/Failed/…`

Depends on 6. One format shared by queue and history, which also settles 5.

## 8. No fixed NCE for Altimeter libraries

- [ ] Suppress the NCE line for Altimeter in the library panel.

`describe_config` in `Routines/DownloadSpecLib/catalog.jl` emits `NCE` from
`nce_params`. The catch: `config.json` does not record the prediction model —
it is known only from `libraries.json`, which is not published yet. So this
needs either the manifest in place, or a rule inferring the model.

## 9. Delay before hover descriptions appear

- [ ] Add an open delay so tooltips do not fire while the pointer crosses them.

`InfoDot.tsx` has no timer; descriptions appear immediately.

## 10. User-defined run name — **done**

- [x] Optional name, generated adjective-noun as the fallback.

Landed in `2b7df871d`. A name already in use gains the next free suffix
(`analysis` → `analysis-2`), checked against persisted history, with the
resolved name previewed under the field.

## Noted while reviewing, not raised

- [ ] The thread picker shows on the DownloadSpecLib page, where it does
      nothing: the transfer is network-bound and the binary ignores it.
