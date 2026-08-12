# Pioneer Console — remaining work

Checklist derived from the meeting plan (`gui_meeting_plan.pdf`) and the design
notes (`gui_design_notes.pdf`). Tick items as they land.

Last updated 2026-08-12.

---

## Next up

- [ ] **C5 — drag a folder onto the browse bar**
      Tauri v2 emits file-drop events per window. Build it once as a shared drop
      target so all three forms get it, with the hover affordance and rejection
      of the wrong kind (a file where a folder is wanted).

- [ ] **D1 — show the library's metadata**
      Decided: read the `.poin`'s own config, not its filename. A filename
      breaks the moment a library is renamed and cannot report anything the
      name does not already say.
      Rust-side read (`paths.rs` already detects `.poin` layout via
      `PION_MARKERS`), then a summary under the field: prediction model, digest
      parameters, modifications, precursor count. Turns the path from a string
      into something checkable before a two-hour run.

- [ ] **F2 — move history and queue to SQLite**
      Needs `tauri-plugin-sql` or `rusqlite`, a schema, and a migration that
      imports whatever is in `localStorage` so existing history survives.
      Do this **before** F3 and F4 — both are much cheaper on top of a real
      store, and doing them first means writing them twice.

- [ ] **F4 — persist the queue across restarts**
      Decided: **record interrupted runs in the history as "interrupted"**, not
      re-queue them. Re-queuing risks starting a long job the moment the app
      opens.
      Their form state must be preserved so the parameters can be recalled and
      re-run without filling the forms out again. `JobSnapshot` already carries
      exactly that for all three commands, and history already restores it — so
      this is mostly a matter of writing pending runs out on exit with a new
      status rather than dropping them.
      Needs an `interrupted` member on `JobStatus`, its own dot colour and
      label, and a re-run affordance on the row.

---

## Small, ready whenever

- [ ] **B3 — "Isotopes (n)"**
      Ambiguous between precursor and fragment isotopes; the parameter is
      fragment isotopes. One line: rename, and add `info` to its `NumSpec` the
      way `Initial NCE` now does.

- [ ] **E3 — threads for BuildSpecLib** *(a measurement, not a feature)*
      Does thread count actually affect a library build? One build at 4 threads
      against one at 15 settles it. The Prosit build measured earlier spent
      300 s of 431 s in Fragment Prediction, so the answer is probably yes —
      which would mean the default should stay high, not drop.

- [ ] **History diff encoding**
      Store only what differs from the defaults. Measured: SearchDIA 528 B →
      278 B (−47%), BuildSpecLib 959 B → 280 B (−71%); the BuildSpecLib saving
      is mostly the four UniProt header regexes, which `presetId: 'uniprot'`
      already implies. `target` also duplicates `results` and can go.
      Independent of F2 — worth doing either way.
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
and configurable, C4 remember the last library, E1 Windows console window,
E2 open the output folder, F1 stable history numbers and no cap.

Also landed, not from either document — the Koina modification picker, the
settings panel behind the gear, and the threads field being clearable.
