Pioneer Search Plan — September 17, 2025

This folder contains ready-to-edit params.json files to run Pioneer SearchDIA across your datasets and Sep 17 libraries. Files are grouped by experiment. Replace the placeholders below, then run the commands in `commands_by_experiment.txt`.

Paths are prefilled for your Windows root:
- ms_data roots: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\{EWZ_MBR | MTAC | OlsenAstralThreeProteome200ng}
- results root: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\Results
- library root: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\SPEC_LIB\\September-17-2025

Notes
- Libraries are the “.poin” directories copied from Sep 17. Do not include a file extension beyond the .poin directory name.
- Defaults are from `assets/example_config/defaultSearchParamsSimplified.json`. Variants override only the needed keys (combine_traces, n_isotopes, huber override, match_between_runs).
- MTAC: adds SeparateTraces for the standard library.
- Proteomes (Yeast, Human, YeastHuman): standard library has extra runs for fragment isotopes (1, 3) and huber delta overrides.
- EWZ_MBR: adds a run with Match-Between-Runs disabled for the standard library.

If any experiment folder names differ (e.g., MTAC is named differently), update the "paths.ms_data" fields in the relevant JSONs.
