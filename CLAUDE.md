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
python scripts/analysis.py   # plots + CTE from output/seq_results.csv (runs at import)
```

Run from the repo root: modules import each other as `scripts.<name>`. There are no tests or
linter. GSAS-II is usually not installed in cloud sessions; logic that doesn't need a real
refinement can be checked by putting an empty `G2script.py` stub on `PYTHONPATH` and faking the
`gpx`/`seqref()` objects.

Inputs: `data/<AFP|WT>/<sample>/*.cbf`, file stems ending in `_<temperature>` (integer K).
`data/` is git-ignored except `data/ice.cif`. Outputs: `projects/<sample>` (.gpx) and `output/`
(timestamped `seq_results_*.xlsx` with Results + Fit_quality sheets, and `seq_results.csv`),
both git-ignored and not auto-created.

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
6. `additional_seq_refine` – sequential again with DisplaceX, background and Mustrain added.
7. `results.extract_data` → `flag_fit_quality` → `save_data`.

Lattice parameters in the sequential fits come from Dij on top of the fixed cell (Cell flag is off);
`seq.get_cell_and_esd` combines them. Variable names are `0:<hist>:D11`, `0:<hist>:Mustrain;a`, etc.

## GSAS-II gotchas

- Refinement recipes are lists of dicts. A repeated key in one dict literal (`'set'`, `'Mustrain'`)
  silently keeps only the last one; this caused real bugs here. One key per step.
- Flags set in a step stay on for later steps until cleared. `'skip': True` applies flags without
  refining.
- Uniaxial Mustrain: `'refine': True` refines both terms; a string refines only that one and clears
  the other; a list raises an error.
- Sequential refinement stops at the first failed histogram; later ones have no entry in
  `seq.data` (`extract_data` handles this).
- GSAS-II source is the reference for what scriptable calls do (`GSASIIscriptable.py`,
  `GSASIIstrMain.py`); check it rather than guessing.

## Physics decisions (see TODO.md for reasoning)

- Samples are live larvae at the tip of a glass capillary: position can change between
  temperatures, and the glass is in the beam.
- LaB6 calibration was collected separately, so absolute a/c carry a distance-related scale error;
  CTE (relative slope) is unaffected if the geometry stayed constant.
- Zero is fixed at 0. Over 6.5–18.5° 2θ it is nearly indistinguishable from DisplaceX, so DisplaceX
  is the only peak-offset parameter (refined per temperature). Don't refine both.
