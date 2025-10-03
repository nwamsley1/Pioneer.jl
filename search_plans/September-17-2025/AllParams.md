Overview of params (summary)

This is a concise summary for quick verification. It shows the three path fields and any non-default overrides per JSON. For full content, open the files under `json/`.

Legend of defaults
- match_between_runs: true
- combine_traces: true
- quant_search.fragment_settings.n_isotopes: 2
- huber_override.override_huber_delta_fit: false

## OlsenAstralThreeProteome200ng (3P OlsenAstral)

- All JSONs
  - ms_data: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\OlsenAstralThreeProteome200ng\\
  - results: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\OlsenAstralThreeProteome200ng\\OlsenAstralThreeProteome200ng_<RunID>
  - variants: Standard; ARAB; EntrapR1/2/3; Standard__Isotopes1; Standard__Isotopes3; Standard__HuberDelta_{100|500|1500|5275|26375|1000000}
  - notes: Isotopes set n_isotopes to 1 or 3; Huber variants set override=true and huber_delta accordingly

## EWZ_MBR (YeastHuman Eclipse)

- All JSONs
  - ms_data: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\EWZ_MBR\\
  - results: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\EWZ_MBR\\EWZ_MBR_<RunID>
  - variants: Standard; ARAB; EntrapR1/2/3; Standard__Isotopes1; Standard__Isotopes3; Standard__HuberDelta_{100|500|1500|5275|26375|1000000}; Standard__MBRoff

## MTAC Three-Proteome (Alternating, Standard)

- All JSONs
  - ms_data: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\MtacThreeProteomeAlternating\\ (or MtacThreeProteomeStandard\\)
  - results: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\<Parent>\\<Parent>_<RunID>
  - variants: Standard; ARAB; EntrapR1/2/3; Standard__Isotopes1; Standard__Isotopes3; Standard__HuberDelta_{100|500|1500|5275|26375|1000000}; Standard__SeparateTraces (combine_traces=false)

## MTAC Yeast (Alternating3M/5M; Standard3M/5M)

- All JSONs
  - ms_data: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\<Parent>\\
  - results: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\<Parent>\\<Parent>_<RunID>
  - variants: Standard; ARAB; EntrapR1/2/3; Standard__Isotopes1; Standard__Isotopes3; Standard__HuberDelta_{100|500|1500|5275|26375|1000000}; Standard__SeparateTraces (combine_traces=false)

## OlsenExplorisThreeProteome500ng (3P Exploris)

- All JSONs
  - ms_data: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\OlsenExplorisThreeProteome500ng\\
  - results: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\OlsenExplorisThreeProteome500ng\\OlsenExplorisThreeProteome500ng_<RunID>
  - variants: Standard; ARAB; EntrapR1/2/3; Standard__Isotopes1/3; Standard__HuberDelta_{...}

## SciexThreeProteomeNswath4 (3P Sciex)

- All JSONs
  - ms_data: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\SciexThreeProteomeNswath4\\
  - results: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\SciexThreeProteomeNswath4\\SciexThreeProteomeNswath4_<RunID>
  - variants: Standard; ARAB; EntrapR1/2/3; Standard__Isotopes1/3; Standard__HuberDelta_{...}

## YEAST_KO (Yeast Olsen)

- All JSONs
  - ms_data: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\YEAST_KO\\
  - results: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\YEAST_KO\\YEAST_KO_<RunID>
  - variants: Standard; ARAB; EntrapR1/2/3

## SingleCell (Human SCP)

- Data and results base
  - ms_data: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\SingleCell\\<Amount>\\<Time>\\
  - results: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\SingleCell\\SingleCell_<Amount>_<Time>_<Variant>
- Libraries
  - Standard: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\SPEC_LIB\\September-17-2025\\altimeter_human_len7o40_ch2o3_mc1_SCP_Sep17-2025.poin
  - ARAB: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\SPEC_LIB\\September-17-2025\\altimeter_human_len7o40_ch2o3_mc1_SCP_Sep17-2025_ARAB.poin
  - EntrapR{1|2|3}: C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\SPEC_LIB\\September-17-2025\\altimeter_human_len7o40_ch2o3_mc1_SCP_Sep17-2025_entrapR{1|2|3}.poin
- Variants by time
  - 24ms (for 10ng, 250pg, 500pg, 1000pg): Standard; ARAB; EntrapR1; EntrapR2; EntrapR3
  - Non-24ms (3ms, 6ms, 12ms): Standard only
  - MBR example (24ms): Standard and Standard__MBRoff (match_between_runs=false)
- Defaults
  - match_between_runs: true (except Standard__MBRoff)
  - combine_traces: true
  - quant_search.fragment_settings.n_isotopes: 2
  - huber_override.override_huber_delta_fit: false
