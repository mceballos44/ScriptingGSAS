# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Automated GSAS-II (scriptable API, `G2script`) sequential Rietveld refinement of ice Ih in frozen
larvae (AFP vs WT samples) over a temperature series, to get a/c lattice parameters, microstrain
and the coefficient of thermal expansion (CTE). Owner: Mauricio Ceballos, Joester Group, Northwestern.

`TODO.md` is the running list of known issues, decisions and their solutions. Read it before
changing refinement logic, and update it (tick items, add new ones with the reasoning) as part of
any change.

## Working rules

- The user wants to understand every change: explain what changed and why, and show the code.
  Commit only when asked; never push to a branch other than the one assigned for the session.
- `data/` (raw `.cbf` images) and `controls/` (calibration, mask, instrument files) hold originals:
  read only. Never edit, move or regenerate them.
- Keep experiment metadata and parameters (sample lists, paths, thresholds, refinement settings)
  separate from code, in config/metadata files. Some still live in `main.py` and `scripts/config.py`;
  moving them out is a TODO item, so don't add new hard-coded parameters.

## Running

Environment is a conda env (`PCScripting`) with `requirements.txt`, plus GSAS-II installed from a
git clone (`pip install ~/G2[gui,useful]`; binaries via `install_gsas_binaries.py`). See
`Instructions.ipynb`.

```
python main.py            # full pipeline for every sample in main.py's `samples` list
python -m scripts.analysis   # plots + CTE from output/seq_results.csv
python -m scripts.stats      # per-sample values + group tests (stats_*.csv)
python -m scripts.report     # stats + output/report.html (d3 from cdnjs)
```

Run from the repo root: modules import each other as `scripts.<name>`. There are no tests or
linter. GSAS-II is usually not installed in cloud sessions; logic that doesn't need a real
refinement can be checked by putting an empty `G2script.py` stub on `PYTHONPATH` and faking the
`gpx`/`seqref()` objects.

Inputs: `data/<AFP|WT>/<sample>/*.cbf`, file stems ending in `_<temperature>` (integer K).
`data/` is git-ignored except `data/ice.cif`. Outputs: `projects/<sample>.gpx` and `output/`
(timestamped `seq_results_*.xlsx` with Results, Fit_quality and Samples sheets, plus
`seq_results.csv` and `sample_info.csv`),
both git-ignored; `config.py` creates the folders on import.

## Pipeline (main.py `full_analysis`, one project per sample)

1. `config.setup` – new project with the ice phase.
2. `image_processing.integrate_images` – import every .cbf sorted by temperature, apply
   `controls/.imctrl` + `.immask`, integrate, set Temperature, load `x.instprm`.
3. `remove_unfrozen` / `remove_orphan_images` – drop patterns without ice peaks (intensity cutoff).
4. `assign_phase_one` + `refinement.initial_refine` – refine **histogram 0** (lowest remaining T)
   alone; it becomes the template for all others.
5. `first_seq_refine` – link phase to all histograms, copy hist 0's background/instrument/limits,
   HAP values and DisplaceX to all, turn on Dij (HStrain), run sequential with `seqCopy`
   (each histogram starts from the previous result).
6. `additional_seq_refine` – sequential again with background and Mustrain added (DisplaceX only
   if `REFINE_DISPLACEMENT_PER_T`).
7. `results.extract_data` → `flag_fit_quality` → `save_data`. `temperature_info` (called between
   steps 2 and 3) records the temperature run and freezing temperature per sample.

`full_analysis` returns `(results_df, sample_info_dict)`; `main.py` skips and reports samples that
raise, so check the "Failed samples" list at the end of a run.
With `config.QUIET_GSAS` on, everything printed during `full_analysis` (GSAS-II and our own
prints) goes to `output/logs/<sample>.log`; only lines matching `GSAS_LOG_KEYWORDS` reach the
terminal. Put user-facing progress messages in `main.py`'s loop, outside the redirect.
Fit-check thresholds live in `config.FIT_CHECKS`; CTE options (freeze scan, excluded scans) in
`config.CTE_SETTINGS`.

Lattice parameters in the sequential fits come from Dij on top of the fixed cell (Cell flag is off);
`seq.get_cell_and_esd` combines them. Variable names are `0:<hist>:D11`, `0:<hist>:Mustrain;a`, etc.

## GSAS-II gotchas

- Refinement recipes are lists of dicts. A repeated key in one dict literal (`'set'`, `'Mustrain'`)
  silently keeps only the last one; this caused real bugs here. One key per step.
- Flags set in a step stay on for later steps until cleared. `'skip': True` applies flags without
  refining.
- Uniaxial Mustrain: `'refine': True` refines both terms; a string refines only that one and clears
  the other; a list raises an error.
- `Rvals['Max shft/sig']` is the total shift from the starting values, not the last cycle; use
  `Rvals['lastShifts']` with `sig` for convergence.
- Histogram `Residuals` (R, wRmin, Durbin-Watson) are only kept from the last single-pattern
  refinement, not per sequential step; compute them from the `data[1]` arrays instead.
- Integrated intensities here are far below counts, so GOF/Rexp are not meaningful in absolute
  terms (GOF ~0.04). Esds are scaled by GOF, so they are still usable.
- Sequential refinement stops at the first failed histogram; later ones have no entry in
  `seq.data` (`extract_data` handles this).
- GSAS-II source is the reference for what scriptable calls do (`GSASIIscriptable.py`,
  `GSASIIstrMain.py`); check it rather than guessing.

## Physics decisions (see TODO.md for reasoning)

- Samples are live larvae at the tip of a glass capillary: position can change between
  temperatures, and the glass is in the beam.
- LaB6 calibration was collected separately, so absolute a/c carry a distance-related scale error;
  CTE (relative slope) is unaffected if the geometry stayed constant.
- Data collection: each sample starts at 260 K or 280 K and is cooled in 1 K steps; every sample
  freezes (around 250 K). Invalid samples are already excluded from the sample list.
- Polarization stays at 0.99 (lab source technician's instruction); don't change it.
- The glass capillary is most of the background. Subtracting the glass image before integration
  made GSAS-II refinements fail, so it's left to the background function for now.
- Zero is fixed at 0. Over 6.5–18.5° 2θ it is nearly indistinguishable from DisplaceX, so DisplaceX
  is the only peak-offset parameter. Don't refine both.
- DisplaceX is refined only in `initial_refine` and then held fixed (`REFINE_DISPLACEMENT_PER_T =
  False`): per-temperature DisplaceX trades off with a and c and doubled their scatter.
- AFP16 and WT20 are excluded (too few frozen scans). The freezing scan is left out of CTE fits,
  and so are individual bad scans in `CTE_SETTINGS['excluded_scans']` (the sample stays in).
- Statistics use one value per larva, never per-temperature points. Start temperature (260/280 K)
  is confounded with collection batch (samples 2-8 vs 15-24).
- Atomic positions and U_iso come from `ice.cif` and are not refined (positions: not determinable
  from this data; U_iso = 0 in the CIF is a TODO).
