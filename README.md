# ScriptingGSAS

Automated GSAS-II analysis of ice in frozen larvae (AFP and WT samples) across a temperature
series. For each sample it integrates the detector images, refines the ice structure at every
temperature, and reports lattice parameters (a, c), microstrain, fit quality and the coefficient of
thermal expansion (CTE).

Mauricio Ceballos, Joester Group, Northwestern University.

## Setup

1. Create and activate the conda environment, then `python -m pip install -r requirements.txt`.
2. Install GSAS-II from a git clone (`git clone --depth 1 https://github.com/AdvancedPhotonSource/GSAS-II.git ~/G2`,
   then `pip install ~/G2[gui,useful]`) and its binaries with `python install_gsas_binaries.py`.
   Details are in `Instructions.ipynb`.
3. Put the images in `data/AFP/<sample>/` or `data/WT/<sample>/`. File names must end in the
   temperature, e.g. `AFP3_250.cbf`.

## Running

```
python main.py               # runs every sample listed in main.py
python -m scripts.analysis   # plots Rwp, a and c vs temperature, prints CTE
python -m scripts.stats      # group comparisons (needs main.py outputs)
python -m scripts.report     # stats + HTML report in output/report.html
```

GSAS-II's printouts for each sample go to `output/logs/<sample>.log`. The terminal shows one
summary line per sample, any GSAS-II line mentioning an error, warning, failure or singular
matrix (prefixed `[GSAS]`), and the summary tables at the end. Set `QUIET_GSAS = False` in
`scripts/config.py` to see everything in the terminal again.

## What each part does

| File | Role |
|---|---|
| `main.py` | Sample list and the per-sample pipeline (`full_analysis`); combines and saves all results |
| `scripts/config.py` | Folder paths and file names; creates a new GSAS-II project with the ice phase |
| `scripts/image_processing.py` | Imports and integrates images, sets temperature, removes unfrozen (no ice peaks) patterns, links the phase |
| `scripts/refinement.py` | Refinement strategy: single-pattern refinement of the first temperature, then sequential refinement of all temperatures |
| `scripts/results.py` | Pulls lattice parameters, strain and fit statistics out of the sequential results, flags questionable fits, writes Excel/CSV |
| `scripts/analysis.py` | Plots and CTE calculation from the saved results |
| `scripts/stats.py` | Per-sample values and statistical tests between groups |
| `scripts/report.py` | Self-contained HTML report: findings, figures and tables (needs internet for fonts and the d3 chart library) |
| `scripts/output_control.py` | Sends GSAS-II printouts to a log file per sample |
| `controls/` | Original calibration and integration files: image controls (`.imctrl`), mask (`.immask`), instrument parameters (`x.instprm`), LaB6 image. Do not edit |
| `data/` | Raw images (not in git) and `ice.cif`. Do not edit |
| `TODO.md` | Known issues, decisions and their status |

## Refinement strategy

1. **First temperature only:** background (5 terms) and sample displacement, then add the unit
   cell and uniaxial microstrain along c. Zero shift stays at 0 (see `TODO.md` 1b).
2. **Copy to all temperatures:** background, instrument parameters, limits, phase settings and
   displacement are copied from the first pattern.
3. **Sequential, pass 1:** refine the Dij strain terms (which give a and c) at each temperature,
   each one starting from the previous temperature's result.
4. **Sequential, pass 2:** also refine background and microstrain. Displacement stays at the value
   from step 1 (refining it at every temperature added noise to a and c; see `TODO.md` 1c).

## Outputs

- `projects/<sample>.gpx`: the GSAS-II project for each sample (open in the GSAS-II GUI to inspect fits).
- `output/seq_results_<date>.xlsx`: sheet *Results* has one row per sample and temperature; sheet
  *Fit_quality* has one row per sample with the flagged temperatures.
  Sheet *Samples* has one row per sample: temperature range and freezing temperature.
- `output/seq_results.csv`: same as *Results*, read by `analysis.py`.
- `output/sample_info.csv`: same as *Samples*.
- `output/logs/<sample>.log`: everything GSAS-II printed for that sample.

### Samples columns

| Column | Meaning |
|---|---|
| `T_start`, `T_end`, `T_step` | Temperature run (cooling, so it starts at `T_start`) |
| `T_freeze` | Warmest scan with ice: the freezing temperature, to within one step |
| `T_last_unfrozen` | The scan just before freezing (no ice yet) |
| `N_unfrozen_below_freeze` | Scans below the freezing point rejected as unfrozen; should be 0, otherwise check those images |

### CTE output (`analysis.calculate_cte`)

Weighted straight-line fit of a or c vs temperature (points with larger error bars count less).
`Red_chi2_fit` is about 1 when the points scatter as much as their error bars say; much larger
means extra scatter or a curved trend (the CTE uncertainty already accounts for it).
Options are set in `CTE_SETTINGS` in `scripts/config.py`: the freezing scan is left out by default,
`excluded_scans` lists individual scans to drop after checking their images, and
`exclude_flagged=True` leaves out every temperature that failed a fit check.

### Statistics (`python -m scripts.stats`)

One value per larva (temperatures of one larva are not independent): CTE of a and c, and a, c
and c/a at a reference temperature (`STATS_SETTINGS['t_ref']`, 248 K) from each sample's
straight-line fit, plus freezing temperature. c/a cancels errors that scale a and c together
(sample-to-detector distance), so it is the most robust structural number.

Tests: freezing temperature 280 vs 260 K start (each group, and both together); AFP vs WT for
every quantity, shuffling labels only within the same start temperature. Each row of
`output/stats_tests.csv` gives group means and medians, the difference with a bootstrap 95%
interval, a permutation p-value, a Mann-Whitney p-value and a Holm-adjusted p-value (corrected
for running all the tests). Per-sample values are in `output/stats_per_sample.csv`.

### Fit-quality columns

| Column | Meaning |
|---|---|
| `Rwp` | Weighted profile R-factor (%): overall misfit between calculated and observed pattern. Compare within a sample, not as an absolute |
| `Rp` | Unweighted profile R-factor (%) |
| `GOF` | Goodness of fit. Assumes intensities in counts; ours are much smaller, so it sits far below 1. Compare between temperatures only |
| `Durbin_Watson` | About 2 when the misfit is random; lower means a systematic misfit (peak shape or background). Finely sampled patterns sit well below 2, so it is flagged only when unusually low for the sample |
| `Converged`, `Last_shift_esd`, `SVD_singular` | Whether the fit settled (final-cycle shift/esd below 1), and whether any parameters were undetermined |
| `Last_shift_param` | The parameter with the largest final-cycle shift/esd |
| `Total_shift_esd` | How far parameters moved from their starting values, in esds (informational, not a convergence test) |
| `Rwp_increase` | Rwp went up compared with the previous temperature |
| `A_trend_z`, `C_trend_z` | How far a or c sits from a straight-line fit vs temperature |
| `DisplaceX` | Sample displacement (fixed after the first refinement by default) |
| `Fit_flags`, `Fit_OK` | Every check that failed, and whether the row passed all of them |
