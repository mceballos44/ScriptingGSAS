# ScriptingGSAS to-do list

Sample setup: live larvae mounted at the end of a glass capillary, lab Ag source (to confirm),
area detector integrated over 300-625° azimuth, 6.5-18.5° 2θ.
LaB6 calibration was collected separately from the sample data.

Status: `[ ]` to do, `[~]` in progress / partly done, `[x]` done

---

## Fit validation

- [x] Add fit metrics to `extract_data` (Rp, Rexp, GOF, χ², Durbin-Watson, convergence, shift/esd, singularities)
- [x] Flag temperatures where Rwp went up from the previous step (`flag_fit_quality`)
- [x] Per-sample summary sheet in the Excel output (`quality_summary`)
- [x] `plot_rwp` in `analysis.py` (circles = Rwp up, x = failed check)
- [ ] Run on real data and pick a sensible `rwp_tol` (default 0 flags every rise)

## 1. Refinement problems that change results

- [x] **1a** `initial_refine`: duplicate dict keys dropped the axial Mustrain setting (step 2) and the
      DisplaceX/Background step (step 5). Rewritten with one key per step.
- [x] **1b** Stop refining Zero on the sample. In 6.5-18.5° 2θ, Zero and DisplaceX shift peaks almost
      identically (cos2θ only goes 0.994 to 0.948), so refining both makes the split arbitrary.
      Decision: Zero fixed at 0, DisplaceX is the only peak-offset parameter.
- [ ] **1c** DisplaceX refined per temperature together with Dij in `additional_seq_refine`.
      Keep for now: a larva on a capillary tip can move between temperatures, so a per-temperature
      offset is physically reasonable. Check that it doesn't fake expansion:
  - [ ] Plot refined DisplaceX vs T per sample (add DisplaceX to `extract_data`)
  - [ ] Compare CTE with DisplaceX fixed vs refined; if they differ, they are correlated
- [ ] **1d** `calculate_cte`: use a weighted fit (`np.polyfit(..., w=1/sigma, cov=True)`) and
      optionally leave out `Fit_OK == False` points
- [ ] **1e** Glass background: subtraction is commented out. The glass capillary sits in the beam,
      so either restore image subtraction or check that 5 background terms are enough
      (low Durbin-Watson = not enough)

## Calibration / sample geometry

- [ ] Sample-to-detector distance: larva at the capillary tip is not where the LaB6 was. This scales
      a and c together and can't be fixed by Zero or DisplaceX. CTE is unaffected if the distance is
      constant; absolute a, c are not.
  - [ ] Compare a, c at one temperature with literature ice Ih (Röttger et al., 1994)
  - [ ] If possible, measure a standard in the same capillary-tip position for future runs
- [ ] Polarization is 0.99 in `.imctrl` and `x.instprm`. If this is a lab Ag source without a
      monochromator, it should be about 0.5 (affects intensities/Rwp, not positions)
- [ ] Larvae can move or dehydrate during a run, and ice in tissue can be spotty (few large crystals).
      Watch for single-temperature Rwp jumps and check those images by eye

## 2. Code that will crash or misbehave

- [ ] `analysis.py` reads `seq_results.xlsx`, which `save_data` never writes (point it at the csv)
- [ ] `analysis.py` runs plots at import; wrap in `if __name__ == "__main__":`
- [ ] `config.py`: create `projects/` and `output/` if missing (`mkdir(exist_ok=True)`)
- [ ] `config.py`: give project files a `.gpx` extension
- [ ] `main.py`: `try/except` around each sample so one failure doesn't lose the whole run
- [ ] `integrate_images`: check the list returned by `Integrate()` before taking `[0]`
- [ ] `requirements.txt`: add `pandas` and `openpyxl`; re-save as UTF-8

## 3. Things that could quietly be wrong

- [ ] `remove_unfrozen`: fixed intensity cutoff depends on exposure; consider a peak-to-median test.
      Check which temperatures it drops per sample
- [ ] `remove_unfrozen`: use `get_sample_name()` so PWDR and IMG names match
      (then `remove_orphan_images` is no longer needed)
- [ ] `LoadProfile` replaces the wavelength with `x.instprm`; confirm it matches `.imctrl` (both 0.56083 now)
- [ ] `copy_displaceX`: copy the list (`list(displace_x)`) instead of sharing one object
- [ ] Size limits in `first_seq_refine` do nothing while Size isn't refined

## Process / project setup

- [ ] Fix GitHub access (Claude GitHub App on `mceballos44/ScriptingGSAS`) and push the local commits
- [x] `CLAUDE.md` with working rules for Claude sessions
- [x] `README.md` documenting what each part does
- [ ] Move parameters out of code into config/metadata files: sample list (`main.py`), paths and
      file names (`config.py`), refinement settings and thresholds, per-sample metadata
      (group, mount, collection date, notes)
- [ ] Enforce read-only `data/` and `controls/` with deny rules in `.claude/settings.json`
- [ ] On the local machine: run `/fewer-permission-prompts` after a few sessions
- [ ] HTML report per run: Rwp and a/c vs T with flags, CTE table, worst fits
- [ ] Once the pipeline is stable: package "run a new sample" as a skill (`skill-creator`)

## 4. Cleanup

- [ ] Remove unused `columns` list in `results.py`
- [ ] Remove `background_list = 0` in `config.py`
- [ ] `ice_cif` defined in four files; keep one
- [ ] Empty `stats_summary` in `analysis.py`
- [ ] Duplicate imports and outdated header comment in `main.py`
