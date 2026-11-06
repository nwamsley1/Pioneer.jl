# Tested Instruments

Instruments whose raw files have been converted and searched with Pioneer, using public data sets and in-house runs. Other instruments or acquisition methods from these vendors may work but have not been tested.

## Thermo Fisher

`.raw` files, converted with **ConvertRAW** (`pioneer convert-raw`).

| Instrument | Mass analyzer (MS2) | Acquisition tested |
|---|---|---|
| Q Exactive | Orbitrap | DIA, including gas-phase fractionation |
| Q Exactive HF | Orbitrap | DIA |
| Q Exactive HF-X | Orbitrap | DIA |
| Orbitrap Exploris 480 | Orbitrap | DIA |
| Orbitrap Eclipse | Orbitrap | DIA |
| Orbitrap Astral | Astral | DIA, including 250 pg single-cell inputs |

## Bruker

`.d` folders, converted with **ConvertRAW** (`pioneer convert-bruker`).

| Instrument | Mass analyzer (MS2) | Acquisition tested |
|---|---|---|
| timsTOF Pro | TIMS-QTOF | diaPASEF |
| timsTOF Pro 2 | TIMS-QTOF | diaPASEF |
| timsTOF SCP | TIMS-QTOF | diaPASEF |
| timsTOF HT | TIMS-QTOF | diaPASEF, Thin-diaPASEF |
| timsTOF Ultra | TIMS-QTOF | diaPASEF |
| timsTOF Ultra 2 | TIMS-QTOF | diaPASEF, including single-cell inputs |
| timsTOF Ultra AIP | TIMS-QTOF | diaPASEF |

Bruker data must be searched with a library built with the **Bruker timsTOF library** option, which adds ion mobility.

## SCIEX

`.wiff` files with their `.wiff.scan` files, converted with **ConvertRAW** (`pioneer convert-sciex`).

| Instrument | Mass analyzer (MS2) | Acquisition tested |
|---|---|---|
| ZenoTOF 7600 | QTOF | SWATH |
| ZenoTOF 8600 | QTOF | ZT Scan DIA |
| TripleTOF 6600 | QTOF | SWATH |
| TripleTOF 5600+ | QTOF | SWATH |
| TripleTOF 5600 | QTOF | SWATH |

`.wiff2` files are not supported; convert them to mzML with msConvert and use **ConvertRAW** with the mzML format.

## Other Instruments

Centroided `.mzML` files from other instruments can be converted with **ConvertRAW** (`pioneer convert-mzml`).
