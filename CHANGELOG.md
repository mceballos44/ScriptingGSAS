# Change log: ScriptingGSAS, 2026-10-08 to 2026-10-09

Everything changed since the original workflow (commit `d5d96b4`, "Initial GSAS-II scripting
workflow"), what the original code did, why each change was made, the decisions taken, and the
results so far. Written as a stopping point for review and discussion.

To see the exact code differences for any file:

```
git diff d5d96b4 -- scripts/refinement.py        # one file
git diff d5d96b4 --stat                          # list of changed files
git show d5d96b4:scripts/refinement.py           # the original version of a file
```

All work is on the branch `claude/cool-davinci-y96ayh`; `main` still has the original code.

---

## 1. Commit timeline

| Commit | What it did |
|---|---|
| `d5d96b4` | **Original workflow** |
| `7d9ba4e` | Fit-quality metrics and Rwp-increase flags in the results |
| `f36dba8` | Rewrote `initial_refine` (duplicate-key bug, no Zero); added `TODO.md` |
| `8771740` | Added `CLAUDE.md` and `README.md` |
| `7d95573`, `506fd62` | TODO updates (statistics plan, GitHub access) |
| `8cbdc45` | Fixed crashes and lost runs (paths, folders, failed samples, requirements) |
| `baf093d` | DisplaceX recorded per temperature; option to hold it fixed |
| `d69e6ad` | Weighted CTE fit |
| `0fe1df7` | Temperature range and freezing temperature recorded per sample |
| `b85ce1c` | GSAS-II printouts sent to a log file per sample |
| `ff8ad5c` | Fixed shift/esd and residual checks after the first real run |
| `e792d60` | DisplaceX fixed after the first refinement; AFP16, WT20 and freezing scans excluded |
| `bd98e77` | Group statistics; six bad scans excluded |
| `71769e0` | HTML report |
| `85ef9d6` | U_iso setting and refinement; ring-spottiness analysis |

---

## 2. The original workflow (`d5d96b4`)

**Pipeline** (`main.py` → `full_analysis` for each of 22 samples):

1. `config.setup`: new project (file `projects/<sample>`, no extension) with the ice phase from
   `data/ice.cif`. Glass background image prepared but commented out.
2. `image_processing.integrate_images`: every `.cbf` sorted by temperature, image controls and
   mask from `controls/`, `sag` and `det2theta` set to 0, integrated, temperature set,
   `x.instprm` loaded.
3. `remove_unfrozen` (99.9th percentile of intensity < 0.80 → unfrozen) and
   `remove_orphan_images`.
4. `assign_phase_one` + `refinement.initial_refine` on histogram 0: six steps (background +
   DisplaceX; cell + Mustrain; clear; Zero once; everything again; clear).
5. `first_seq_refine`: link phase to all histograms, copy histogram 0's background, instrument
   parameters, limits, HAP values and DisplaceX; refine Dij (HStrain); sequential with `seqCopy`,
   10 cycles, Size limits 0.001-5.
6. `additional_seq_refine`: sequential again with DisplaceX, background and Mustrain.
7. `results.extract_data` → `combine_data` → `save_data` (timestamped `.xlsx` + `seq_results.csv`).
   Columns: sample, T, a, c, D11, D33, Mustrain (a, i) with esds, Rwp.

`scripts/analysis.py` plotted a and c vs T and calculated CTE with an unweighted `linregress`.

### Problems found in the original code

| Where | Problem | Effect | Status |
|---|---|---|---|
| `initial_refine` step 2 | `'Mustrain'` key written twice in one dict; Python keeps only the last | The axial / `[0,0,1]` setting was silently dropped; only the equatorial term refined | Fixed (`f36dba8`) |
| `initial_refine` step 5 | `'set'` key written twice | DisplaceX + background refinement in that step never happened | Fixed (`f36dba8`) |
| `initial_refine` step 4 | Zero refined on the sample | Zero and DisplaceX are nearly identical over 6.5-18.5° 2θ; the split was arbitrary | Zero now fixed at 0 (`f36dba8`) |
| `analysis.py` | Read `output/seq_results.xlsx`, which `save_data` never wrote | Plots and CTE crashed | Fixed (`8cbdc45`) |
| `analysis.py` | Plots ran at import | Importing the module ran everything | Fixed (`8cbdc45`) |
| `config.setup` | `projects/` and `output/` not created; project file had no `.gpx` | Fresh checkout failed; GUI couldn't open projects | Fixed (`8cbdc45`) |
| `main.py` | One failed sample stopped the whole run | All results lost | Fixed (`8cbdc45`) |
| `integrate_images` | `Integrate()[0]` before checking the result | IndexError instead of a clear message | Fixed (`8cbdc45`) |
| `extract_data` | `seq.get_Variable(...)[0]` on missing variables; histograms the sequential refinement never reached | Crash | Fixed (`7d9ba4e`) |
| `requirements.txt` | Missing `pandas`, `openpyxl`; saved as UTF-16 | Install incomplete; git shows it as binary | Fixed (`8cbdc45`) |
| `calculate_cte` | Unweighted fit, ignored esds | Less precise, overconfident CTE errors | Fixed (`d69e6ad`) |
| `remove_unfrozen` | Fixed intensity cutoff; name matching only works because `remove_orphan_images` cleans up | Fragile | Open (TODO §3) |
| `ice.cif` | U_iso = 0 for every atom | Intensities fall off too slowly with angle | Addressed (`85ef9d6`), to be tested |

---

## 3. Changes by file

### `scripts/refinement.py`
- **`initial_refine`** rewritten: three steps, one `set` per step.
  1. background (5 terms) + DisplaceX;
  2. + cell + uniaxial Mustrain along [001] (`'refine': True` = both terms) + O U_iso (if
     `UISO_SETTINGS['mode'] == 'refine'`);
  3. clear all flags without refining (`'skip': True`).
  Zero is no longer refined.
- **`additional_seq_refine`** takes `refine_displacement` (default from
  `config.REFINE_DISPLACEMENT_PER_T = False`): DisplaceX is no longer refined per temperature.

### `scripts/results.py`
- `extract_data`: safe variable lookup; rows for histograms the sequential refinement didn't
  reach (`Refined = False`); new columns `Seq_index`, `Histogram`, `DisplaceX` (+ esd), `Uiso_O`.
- Fit metrics per temperature: Rwp, Rp, GOF, reduced χ², Durbin-Watson, Nobs, Nvars, Converged,
  DelChi2, `Last_shift_esd` + `Last_shift_param`, `Total_shift_esd`, SVD singularities, aborted,
  GSAS-II message. Rp and Durbin-Watson are computed from the pattern arrays (GSAS-II only kept them
  for the first histogram) and checked against GSAS-II's own Rwp.
- `flag_fit_quality`: Rwp rise from the previous temperature, Rwp / GOF outliers, low
  Durbin-Watson for the sample, not converged, last-cycle shift/esd, singular parameters, a or c
  off its straight-line trend. Gives `Fit_flags` and `Fit_OK`.
- `quality_summary` (one row per sample) and `temperature_info` (temperature range, freezing
  bracket, frozen-looking frames rejected below freezing).
- `save_data` writes three sheets (Results, Fit_quality, Samples) plus `seq_results.csv` and
  `sample_info.csv`.

### `scripts/config.py`
Now holds the settings, so they're out of the code:
`QUIET_GSAS`, `LOG_DIR`, `GSAS_LOG_KEYWORDS`, `REFINE_DISPLACEMENT_PER_T`, `UISO_SETTINGS`,
`CTE_SETTINGS` (freezing scan, excluded scans), `FIT_CHECKS` (thresholds), `SPOTTINESS_SETTINGS`,
`STATS_SETTINGS`. Creates `projects/` and `output/`; project files get `.gpx`; `set_uiso` applies
the U_iso starting values.

### `main.py`
- AFP16 and WT20 commented out of the sample list (too few frozen scans).
- Each sample runs inside a log redirect and a `try/except`: GSAS-II output goes to
  `output/logs/<sample>.log`, failures are reported and skipped.
- One summary line per sample (range, freezing temperature, scans refined, Rwp range).
- Fit flags, sample info and summaries are saved and printed at the end.

### `scripts/analysis.py`
- Reads `output/seq_results.csv`; plots only when run directly (`python -m scripts.analysis`).
- `calculate_cte`: weighted `np.polyfit` with `cov=True`; skips unrefined scans, the freezing scan
  and excluded scans; optional `exclude_flagged`; reports points used, `Red_chi2_fit`, and the
  fitted value at a reference temperature.
- `plot_rwp` (Rwp vs T with flags) and a DisplaceX plot.

### `scripts/image_processing.py`
- Only the `Integrate()` result check changed.

### New files
| File | Purpose |
|---|---|
| `scripts/output_control.py` | Sends printed output to a log file, passing lines with warning keywords to the terminal |
| `scripts/stats.py` | One value per larva; permutation tests, bootstrap CIs, Mann-Whitney, Holm |
| `scripts/report.py` | Runs the stats and writes `output/report.html` |
| `scripts/spottiness.py` | Ring spottiness from the raw images (currently fails, see §6) |
| `TODO.md` | Running list of issues, decisions and status |
| `README.md` | What each part does, how to run, what the outputs mean |
| `CLAUDE.md` | Rules and context for Claude sessions in this repo |
| `CHANGELOG.md` | This file |

### `requirements.txt`
Re-saved as UTF-8; `pandas` and `openpyxl` added without versions (re-run `pip freeze` locally to
pin them).

---

## 4. Decisions and the evidence behind them

| Decision | Why |
|---|---|
| Zero fixed at 0; DisplaceX is the only peak offset | Over 6.5-18.5° 2θ, cos 2θ changes only 0.994 → 0.948, so Zero and DisplaceX shift peaks almost identically |
| DisplaceX refined once, then held (`REFINE_DISPLACEMENT_PER_T = False`) | Comparison runs 15:53 (refined) vs 16:03 (fixed): per-temperature DisplaceX scattered less than its esd, showed no drift, and was anti-correlated with the a/c scatter (r = -0.6 to -0.97). Fixing it halved the scatter of a and c and the CTE uncertainty, and removed most of an apparent AFP-WT difference in CTE of a (7 → 2 ×10⁻⁶ K⁻¹) |
| Polarization stays 0.99 | Lab source technician's instruction |
| Glass background left to the background function | Subtracting the glass image before integration made refinements fail |
| AFP16, WT20 excluded | 4 and 2 frozen scans |
| Freezing scan excluded from fits | Partly liquid; AFP5 250 K image is noisy with incomplete rings |
| Scans excluded (sample kept): AFP6 243, AFP17 246, WT23 246, WT3 241/244/245 | Rwp 17-25% vs ~8%, the same with DisplaceX fixed or refined, so the patterns themselves are off |
| Fit-check thresholds: `rwp_tol` 0.05, `max_shift_esd` 1.0 | First run: Rwp noise 0-3% between neighbours; GSAS-II stops with last shifts ~0.7 esd (mostly Mustrain, not a/c) |
| Statistics on one value per larva | Temperatures of one larva aren't independent |
| Compare a, c, c/a at 248 K | Inside every larva's range except WT21 (extrapolated 2 K) |
| O and H positions not refined | d ≥ 1.75 Å (~a dozen reflections), H nearly invisible, spotty intensities; positions don't affect a, c |
| Crystallite size not refined | Size broadening (1/cos θ) is flat over this range, so it can't be separated from constant-width terms; micron grains don't broaden peaks anyway |

---

## 5. Results so far (2026-10-09 16:27 run, before the U_iso change)

| Comparison | Difference [95% CI] | p | Holm |
|---|---|---|---|
| Freezing T, start 280 vs 260 K (both groups) | +2.4 K [1.0, 3.7] | 0.005 | 0.04 |
| Freezing T, AFP vs WT | AFP 1.5 K colder [0.3, 2.7] | 0.06 | 0.38 |
| CTE of a, AFP vs WT | +1.1 [-1.3, 3.4] ×10⁻⁶ K⁻¹ | 0.42 | 1 |
| CTE of c, AFP vs WT | -0.1 [-4.0, 3.5] ×10⁻⁶ K⁻¹ | 0.94 | 1 |
| c/a at 248 K, AFP vs WT | -0.00004 [-0.00015, 0.00006] | 0.47 | 1 |

- CTE of a ≈ 49 and of c ≈ 46 ×10⁻⁶ K⁻¹ in both groups.
- Start temperature coincides with collection batch (larvae 2-8 at 260 K, 15-24 at 280 K), so a
  batch difference would look the same.
- Report for this run: https://claude.ai/artifact/X2JUFAnPwcD12jPoynuNVa (private until shared).

---

## 6. Where things stopped

**Spottiness analysis fails** at the mask step (`ValueError: operands could not be broadcast
together with shapes (301453,) (237169,)` in `GSASIIimage.MakeMaskMap`). The detector is 487 × 619
pixels (301453 = 487 × 619; 237169 = 487 × 487). `MakeMaskMap` builds its polygon mask from
`controls['size']` in GSAS-II's (x, y) order, while the image array is in (row, column) order, so
for a non-square detector the two disagree. Proposed fix: build the polygon/point mask directly
from the pixel coordinates in mm (as `Make2ThetaAzimuthMap` does) instead of calling
`MakeMaskMap`, and test it on one image. Not changed yet.

**Waiting on a run:**
- U_iso: compare `UISO_SETTINGS['mode']` `'refine'` vs `'cif'` (Rwp, a, c, CTE).

**Still open** (details in `TODO.md`): Ag Kα₁/Kα₂ doublet (LaB6 recalibration by hand), preferred
orientation test, empirical peak width, absolute-scale check against literature ice Ih,
`remove_unfrozen` robustness, convergence of Mustrain, moving the sample list into config,
read-only protection for `data/` and `controls/`, code cleanup.

---

## 7. Questions to discuss

1. Start temperature vs batch: were larvae 2-8 and 15-24 collected on different days or from
   different cohorts? Is there a way to separate the two?
2. The AFP vs WT freezing difference (1.5 K, p = 0.06): worth more samples?
3. No difference in lattice or CTE between groups: is that the expected outcome, or does it
   change the focus toward grain size (spottiness)?
4. Should excluded scans be decided by a rule (e.g. Rwp > 2× the sample median) rather than by
   inspection, for reproducibility?
5. Is the LaB6 recalibration with the Kα doublet worth doing before the width analysis?
