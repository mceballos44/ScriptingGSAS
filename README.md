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
4. Create the `projects/` and `output/` folders (not created automatically yet).

## Running

```
python main.py               # runs every sample listed in main.py
python scripts/analysis.py   # plots Rwp, a and c vs temperature, prints CTE
```

## What each part does

| File | Role |
|---|---|
| `main.py` | Sample list and the per-sample pipeline (`full_analysis`); combines and saves all results |
| `scripts/config.py` | Folder paths and file names; creates a new GSAS-II project with the ice phase |
| `scripts/image_processing.py` | Imports and integrates images, sets temperature, removes unfrozen (no ice peaks) patterns, links the phase |
| `scripts/refinement.py` | Refinement strategy: single-pattern refinement of the first temperature, then sequential refinement of all temperatures |
| `scripts/results.py` | Pulls lattice parameters, strain and fit statistics out of the sequential results, flags questionable fits, writes Excel/CSV |
| `scripts/analysis.py` | Plots and CTE calculation from the saved results |
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
4. **Sequential, pass 2:** also refine displacement, background and microstrain.

## Outputs

- `projects/<sample>`: the GSAS-II project for each sample (open in the GSAS-II GUI to inspect fits).
- `output/seq_results_<date>.xlsx`: sheet *Results* has one row per sample and temperature; sheet
  *Fit_quality* has one row per sample with the flagged temperatures.
- `output/seq_results.csv`: same as *Results*, read by `analysis.py`.

### Fit-quality columns

| Column | Meaning |
|---|---|
| `Rwp` | Weighted profile R-factor (%): overall misfit between calculated and observed pattern. Compare within a sample, not as an absolute |
| `Rexp`, `GOF` | Best Rwp the counting noise allows, and Rwp/Rexp (ideally near 1) |
| `Durbin_Watson` | About 2 when the misfit is random; below 1 means a systematic misfit (peak shape or background) |
| `Converged`, `Max_shift_esd`, `SVD_singular` | Whether the fit settled, and whether any parameters were undetermined |
| `Rwp_increase` | Rwp went up compared with the previous temperature |
| `A_trend_z`, `C_trend_z` | How far a or c sits from a straight-line fit vs temperature |
| `Fit_flags`, `Fit_OK` | Every check that failed, and whether the row passed all of them |
