# Thermo .raw reader (pure Julia, no Thermo DLL)

Reverse-engineered from data files only, checked scan by scan against `.arrow` files written by PioneerConverter
(Thermo RawFileReader). Reader: `ThermoRaw.jl` (layout documented at the top). Check: `verify_reader.jl <raw> <arrow>`.
Prototype only, not wired into Pioneer. `probes/` holds the exploration scripts in the order the format was worked out;
they hard-code local copies under `~/BrukerTims/rawformat/` (sources: RIS `NTW/PioneerTutorial/RawData`,
`Raw files/MTAC/Astral Alternating Windows`, and `~/Projects/pioneer-gui-test/rawtest`).

## Verified (format version 66), zero mismatches on every field

| File | Instrument | Scans | Packet types |
|---|---|---|---|
| Olsen Exploris E5H50Y45 DIA_1 | Orbitrap Exploris 480 | 62,196 | 21 (Orbitrap profile + centroids) |
| HELA_MOUSESIL_METHTEST_07312024_06 | Orbitrap Eclipse (ion-trap + Orbitrap) | 122,809 | 18 (ion trap), 21 |
| MTAC Yeast Alternating-v2 3min Rep1 | Orbitrap Astral | 38,020 | 21 (MS1), 20 (Astral MS2, centroids only) |

Fields matched exactly: centroid m/z and intensity, RT, TIC, base peak m/z/intensity, scan low/high m/z, packet type,
MS order, isolation centre, isolation width, collision energy. Opening takes 0.2-2 s (self-pointer scan); reading all
scans about 1 s.

## Behaviour that mirrors the converter
- Peaks flagged 0x10 in the per-peak descriptor word (reference / lock-mass) are dropped.
- Ion-trap packets (type 18) produce no peaks, as with RawFileReader's centroid stream.
- MS1 records also carry a "reaction" (centre/width of the scan range); the converter reports missing for MS1.

## Not done yet
- `collisionEnergyEvField` (likely in the trailer key/value records at RunHeader+7456) and the scan filter string.
- Locate the controller list by parsing the variable-length file info, instead of scanning for the RunHeader
  self-pointer (works, costs ~1.5 s per GB).
- Older format versions (< 66), profile-only scans, and writing Pioneer's .arrow directly.
- Licensing review before shipping in Pioneer (work done from data files only, never the DLL).
