# Instrument test data

Where the raw files behind the [Tested Instruments](../docs/src/user_guide/instruments.md) page live, for re-testing
conversion and search and for adding regression sets. Results and status are as of October 2026.

RIS paths are relative to `/storage3/fs1/d.goldfarb/Active/`; `NTW` is `RIS_Goldfarb_Lab/NTW`.
`NTW/INSTRUMENT_TESTDATA_2026-10/` holds raw copies of files that had no other RIS copy, plus `provenance.tsv`
(accession, source URL, paper, RIS copy and result for every file) and `md5sums.txt`.

## Thermo (`.raw` → PioneerConverter → `.arrow`)

| Instrument | MS2 analyzer | Accession | Test file | RIS copy | Status |
|---|---|---|---|---|---|
| Q Exactive | Orbitrap | PASS00589 | `B_D140314_SGSDSsample6_R01_MHRM_T0.raw` | `NTW/INSTRUMENT_TESTDATA_2026-10/thermo/QE_PASS00589.raw` | searches |
| Q Exactive (GPF) | Orbitrap | PXD055536 | `GPF6_HM_POOL_DIA_6.raw` | `NTW/INSTRUMENT_TESTDATA_2026-10/thermo/QE_GPF_PXD055536.raw` | searches |
| Q Exactive HF | Orbitrap | PXD005573 | `Fig2HeLa-0-5h_MHRM_R01_T0.raw` | `NTW/INSTRUMENT_TESTDATA_2026-10/thermo/QEHF_PXD005573.raw` | searches |
| Q Exactive HF-X | Orbitrap | PXD028735 | `LFQ_Orbitrap_AIF_Ecoli_02.raw` | `NTW/INSTRUMENT_TESTDATA_2026-10/thermo/QEHFX_15k_PXD028735.raw` | searches |
| Orbitrap Exploris 480 | Orbitrap | PXD073390 | `E23501_MS24-055_970-2_MAM-GSC-9-3_12.raw` | `NTW/INSTRUMENT_TESTDATA_2026-10/thermo/Exploris480_PXD073390.raw` | searches |
| Orbitrap Exploris | Orbitrap | PXD046444 | Olsen three-species mix | `Automation/Pioneer/data/Olsen_Exploris_3P` (arrows) | regression set |
| Orbitrap Astral | Astral | PXD046444 | Olsen three-species mix | `Automation/Pioneer/data/Olsen_Astral_3P` (arrows) | regression set |
| Orbitrap Astral | Astral | in-house | MTAC three-proteome standard | `Automation/Pioneer/data/MTAC_3P_Standard` (arrows) | regression set |
| Orbitrap Astral | Astral | not recorded | 250 pg, 3/6/12/24 ms IT, FAIMS | `Automation/Pioneer/data/SCP_Astral_250pg_*` (arrows) | regression set |
| Orbitrap Astral | Astral | PXD070049 | `LFQ_Astral_DIA_5min_250pg_Ecoli_01.raw` | `NTW/INSTRUMENT_TESTDATA_2026-10/thermo/Astral_PXD070049.raw` | searches |
| Orbitrap Eclipse | Orbitrap | PXD023334 | `20201024_LRH_plasma_OT_1-7_2.raw` | `NTW/INSTRUMENT_TESTDATA_2026-10/heil_PXD023334/` | searches |
| Orbitrap Eclipse | Orbitrap | PXD038508 | `20210906_LRH_MMCC_Static_K_1.raw` (raw metadata says Eclipse) | `NTW/INSTRUMENT_TESTDATA_2026-10/thermo/Lumos_15k_PXD038508_K_1.raw` | searches |
| Orbitrap Eclipse | Orbitrap, 1 ng | MSV000089718 | `20220220_EV_Evo_Whisper100_40SPD_HeLa_1ng_1CV45_OTDIA_15k_s01.raw` | `NTW/INSTRUMENT_TESTDATA_2026-10/thermo/Eclipse_OTDIA15k_MSV000089718.raw` | weak (1 ng) |
| Orbitrap Eclipse | ion trap | PXD023334 | `20201024_LRH_plasma_IT_1-7_2.raw` | `NTW/INSTRUMENT_TESTDATA_2026-10/heil_PXD023334/` | not supported yet |
| Orbitrap Eclipse | ion trap Turbo, 1 ng | MSV000089718 | `20220221_EV_Evo_Whisper100_40SPD_HeLa_1ng_1CV45_LITDIA_Turbo_8ms_s01.raw` | `NTW/INSTRUMENT_TESTDATA_2026-10/thermo/Eclipse_LITDIA_Turbo_MSV000089718.raw` | not supported yet |
| Orbitrap Eclipse | ion trap Rapid, 1 ng | MSV000090688 | `20220708_HB_ES_Evo_Whisper100_40SPD_LIT-DIA_40w_10mz_Rapid_auto_1CV_45_1ng_1.raw` | `NTW/INSTRUMENT_TESTDATA_2026-10/thermo/Eclipse_LITDIA_Rapid_MSV000090688.raw` | unsearchable (accepted) |
| Stellar | ion trap | MSV000094904 / PXD052847 | see below | `NTW/INSTRUMENT_TESTDATA_2026-10/stellar_MSV000094904/` | not supported yet |

## Bruker (`.d` → `convertBruker` → `.tdfs`)

| Instrument | Acquisition | Accession | Test file | RIS copy | Status |
|---|---|---|---|---|---|
| timsTOF Pro | diaPASEF | PXD040205 | `20220715_1000ngK562_ZI_500ul60C3cm_3min_P1-D1_1_9772.d` | `BrukerTims_Datasets/PXD040205` | searches |
| timsTOF Pro | diaPASEF | PXD053462 | `20220407_KK_HeLa_Ecoli-10perc_10ng_R2_F5_1_6903.d` | `BrukerTims_Datasets/PXD053462` | searches |
| timsTOF Pro 2 | diaPASEF | PXD029836 | `20211101_PRO2_LS_12_MA_HeLaSCS_10.0_ngHS_GE4_1_1416.d` | `BrukerTims_Datasets/PXD029836` | searches |
| timsTOF SCP | diaPASEF | PXD062685 | `ttSCP_diaPASEF_Condition_B_Sample_Alpha_01_11496.d` | `BrukerTims_Datasets/PXD062685` | searches |
| timsTOF SCP | diaPASEF, low input | PXD044986 | `EGI_C13_DIA_IO_L15_Slot2-41_1_1071.d` | `NTW/INSTRUMENT_TESTDATA_2026-10/bruker/timsTOF_SCP_PXD044986.d.tar.gz` | searches |
| timsTOF HT | Thin-diaPASEF | PXD061613 / JPST003679 | `HEK100ng_Thin-diaPASEF-1.d`, `H65Y15E20_200ng_Thin-diaPASEF-3.d` | `BrukerTims_Datasets/PXD061613` | searches |
| timsTOF HT | slice-PASEF | PXD061613 / JPST003679 | `HEK100ng_Ramp100ms_075-120_slicePASEF-1.d` | `NTW/SLICEPASEF_HANDOFF_2026-10-19/raw` | not supported yet |
| timsTOF HT | diaPASEF | PXD078816 | `6743_HT3_DIA_RSLC_2col-108min_35_1_BC3_1_11347.d` | `NTW/INSTRUMENT_TESTDATA_2026-10/bruker/timsTOF_HT_PXD078816.d.zip` | searches |
| timsTOF Ultra | diaPASEF | PXD071819 | `DIA_SampleA_1ng_S1-A2_1_12107.d` | `BrukerTims_Datasets/PXD071819` | searches |
| timsTOF Ultra | diaPASEF | PXD070049 | `LFQ_Ultra_diaPASEF_5min_250pg_Condition_B_REP6.d` | `NTW/INSTRUMENT_TESTDATA_2026-10/bruker/timsTOF_Ultra_PXD070049.d.zip` | searches |
| timsTOF Ultra 2 | diaPASEF | PXD070049 | `LFQ_Ultra2_diaPASEF_15min_50ng_Condition_A_REP1.d` | `NTW/PXD070049_Bruker` | searches |
| timsTOF Ultra 2 | diaPASEF | PXD075483 | `250314_HeLa_100ng_90min_DIA_01_S2-A1_1_507.d` | `BrukerTims_Datasets/PXD075483` | searches |
| timsTOF Ultra 2 | diaPASEF, single cell | PXD079517 | `H358_control_diaPASEF_80SPD_S2-H3_1_476.d` | `BrukerTims_Datasets/PXD079517` | searches |
| timsTOF Ultra AIP | diaPASEF | PXD069877 | `Buttercup_20250410_MKL_P03E15_120SPD_Tonsil_8k_3x...d` | `BrukerTims_Datasets/PXD069877` | searches |

## SCIEX (`.wiff` + `.wiff.scan` → `convertSciex` → `.scxs`)

| Instrument | Acquisition | Accession | Test file | RIS copy | Status |
|---|---|---|---|---|---|
| ZenoTOF 7600 | nSWATH | not recorded | `BenchSample{A,B}_nswath4_Rep{1,2,3}` | `Automation/Pioneer/data/SCIEX_3P` (arrows) | regression set |
| ZenoTOF 8600 | ZT Scan DIA | PXD070049 | `LFQ_ZenoTOF8600_ZTScanDIA_5Da_30min_Nano_50ng_Condition_A_REP1.wiff` | `Automation/Pioneer/data/SCIEX_ZT`; `NTW/SciexZT_ZTScanDIA` | regression set |
| TripleTOF 6600 | SWATH | PXD014690 | `20181206_k562_2ug_20min_3.wiff` | `NTW/INSTRUMENT_TESTDATA_2026-10/sciex/TTOF6600_K562_PXD014690.wiff(.scan)` | searches |
| TripleTOF 5600+ | SWATH | PXD010912 | `HLS9_DIA33.wiff` | `NTW/INSTRUMENT_TESTDATA_2026-10/sciex/TTOF5600plus_PXD010912.wiff(.scan)` | searches |
| TripleTOF 5600 | SWATH | PXD028735 | `LFQ_TTOF5600_SWATH_Ecoli_01.wiff` | `NTW/INSTRUMENT_TESTDATA_2026-10/sciex/TTOF5600_PXD028735.wiff(.scan)` | searches |
| TripleTOF 5600 | SWATH | PXD001587 | `LongSwath_UPS1_40fm_Ecoli_1ug-rep3.wiff` | `NTW/INSTRUMENT_TESTDATA_2026-10/sciex/TTOF5600_PXD001587.wiff(.scan)` | weak (sample) |
| TripleTOF 6600 | SWATH | PXD028735 | `LFQ_TTOF6600_SWATH_Ecoli_03.wiff` | `NTW/INSTRUMENT_TESTDATA_2026-10/sciex/TTOF6600_PXD028735.wiff` | unreadable: deposit has no `.wiff.scan` |

## Thermo Stellar

Searle lab GPF-DIA data, mouse immune cells (Nat Commun 2025; MassIVE MSV000094904, PRIDE PXD052847): 18 raw files
(IL2A, IL15_PRTC and ILPool_PRTC, each in 6 gas-phase fractions of 100 m/z), 2 Th windows, 60-min gradient.

- **RIS copy:** `NTW/INSTRUMENT_TESTDATA_2026-10/stellar_MSV000094904/` holds
  `20240408_P5_Neo_ES75150_500ng_{IL2A_GPFDIA_60min_02,IL15_PRTC_GPFDIA_60min_04,ILPool_PRTC_GPFDIA_60min_04}.raw`.
- **Corrupt public copies:** the MassIVE copies of IL2A_03, IL2A_06, IL15_06, ILPool_03 and ILPool_05 are corrupt
  (their SHA-1 differs from the source hash msConvert recorded in the authors' mzML). Download those from Panorama
  Public instead.
- **Conversion:** PioneerConverter v1.2.0 writes zero peaks for ion-trap (ITMS) scans. Converting Stellar data needs
  PioneerConverter PR #14 ("Support Thermo Stellar (ion-trap) raw files").
- **Library:** the prebuilt Altimeter *Mus musculus* library (`ListSpecLibs()` / `pioneer download --list`).
- **October 2026 result:** 9,900–13,800 precursors and 3,800–4,200 protein groups per file, with PR #14 and the
  unit-resolution tolerance changes. Re-test on the current develop before adding Stellar to the Tested Instruments
  page or a regression set.
