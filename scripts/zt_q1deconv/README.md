# Library-free Q1 deconvolution of ZT scanning-DIA (prototype, 2026-09-18)

Converts a Sciex ZT scanning-quad Arrow file into a conventional 1-Da DIA file that Pioneer
searches with `scanning_quad: false`. Per cycle: bin MS2 peaks on an m/z grid, Gaussian-smooth
across grid rows, NNLS-deconvolve each row along Q1 with the measured triangular transmission
kernel (half-base h/S = 6.49 bins, k = 6 edge slices), optionally smooth the solution along Q1
with a 3-tap kernel, emit one spectrum per slice at the intensity-weighted raw m/z.

    JULIA_NUM_THREADS=12 julia --project=. zt_q1_deconvolve.jl <in.arrow> <out.arrow> <rt_lo> <rt_hi> [raw_window.arrow]

Env: `MODE` nnls|centroid|trace|tracecen, `GRID_MDA` (10), `SIGMA_BINS` (1), `H_BINS`, `K_BINS`,
`MIN_INT` (20), `POST_Q1_H` (3-tap [1-1/H, 1, 1-1/H]; 3 is best so far), `POST_Q1_TRI`,
`CD_OMEGA` (1.6), `TRACE_PPM`, `TRACE_MINLEN`.

## Results, A_REP1 RT 6.5-7.5 (56 cycles), Pioneer precursors at 1% FDR
DIA-NN (whole-file first pass, apex RT in window): 11,704. Pioneer ZT path on the raw window: 10,449.

| converter | grid | σ | post-Q1 | precursors | note |
|---|---|---|---|---:|---|
| grid NNLS | 5 mDa | 1 | – | 7,608 | first version |
| grid centroid | 5 | 1 | – | 5,781 | local maxima; loses shoulders |
| trace NNLS (no grid) | – | – | – | 6,971 | links peaks bin-to-bin, 10 ppm |
| trace centroid | – | – | – | 6,123 | |
| grid NNLS | 5 | 2 | – | 8,564 | |
| grid NNLS | 10 | 1 | – | 8,686 | |
| grid NNLS | 10 | 2 | – | 7,830 | |
| grid NNLS | 20 | 1 | – | 7,481 | |
| grid NNLS | 5 | 1 | ½ 1 ½ | 7,300 | |
| grid NNLS | 10 | 1 | ½ 1 ½ | 9,596 | |
| grid NNLS | 10 | 1 | ⅔ 1 ⅔ | **9,704** | 93% of ZT path, 83% of DIA-NN |
| grid NNLS | 10 | 1 | ⅚ 1 ⅚ | 9,016 | too flat: three equal candidates per precursor |

Square quad model instead of the Razo fit on the converted file: no gain (7,521 / 7,648 vs 7,608).
The Razo fit on a converted file is a triangle of half-base ~1 slice, which is what discrete
slices produce for a precursor between two slice centres.

## Findings
- Fragment m/z along a streak is correlated drift, not iid jitter: lag-1 autocorrelation 0.38
  (median), 15 ppm span over 13 bins, 88% of streaks cross more than one 5 mDa row (`q1_mz_drift.jl`).
- Q1 smoothing before NNLS is a no-op when the kernel is made composite; with the original kernel
  it leaks. m/z smoothing is the step that matters (`q1_smooth_panels.jl`).
- Residual artefact of the row-wise grid: drifting streaks become diagonal chains of small
  spurious slices. Coarser rows / post-Q1 smoothing reduce the damage.
- Converter cost after optimisation: ~6 s/cycle on 8 threads (grid 0.9, solve 4.9, emit 0.1),
  ~15 min per 5-min file on 12 threads; was ~2 h. The solve is coordinate descent with ω=1.6;
  an exact Lawson-Hanson per segment is more accurate but needs an allocation-free implementation.
- Search of a converted file: main search 3-6 s vs 46 s for the ZT path on the same window.
