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
- [x] Run on real data and pick a sensible `rwp_tol`: 0.05. First run (2026-10-09): rises of 0-3%
      are noise, real problems were +80-200%. Thresholds now in `FIT_CHECKS` (`config.py`)
- [x] Fix: "large shift/esd" flagged 244/246 rows. GSAS-II's `Max shft/sig` is the total change
      from the starting values, not the last cycle. Now `Last_shift_esd` from `Rvals['lastShifts']`
- [x] Fix: Rp, Rexp and Durbin-Watson were only stored for the first histogram (stale values from
      `initial_refine`). Rp and DW are now computed per temperature from the pattern arrays
      (checked against GSAS-II's Rwp); Rexp dropped. DW flagged relative to the sample, since
      every pattern sits around 0.4-1.1
- [ ] GOF ~0.04 and Rexp = 100%: integrated intensities are far below counts (the unfrozen cutoff is
      0.8), so GSAS-II's 1/y weights don't match the real noise. Rwp and esds (scaled by GOF) are
      fine; GOF/Rexp are not interpretable in absolute terms. Find where the scaling comes from
      (image controls / detector gain) if absolute GOF matters
- [x] Bad fits in the first run (Rwp jumps from ~8% to 17-25%), check these patterns and images.
      Decision: leave these scans out of the analysis (not the samples); listed in
      `CTE_SETTINGS['excluded_scans']`.
      They are identical with DisplaceX fixed, so the patterns themselves are off, not the fit.
      AFP5 250 checked: noisy image, incomplete rings (freezing artifact). Once checked, list bad
      scans in `CTE_SETTINGS['excluded_scans']` (`config.py`):
      AFP16 247, AFP17 246, AFP5 250, AFP6 243, AFP8 250, WT23 246, WT3 241/244/245.
      AFP5 250 and AFP8 250 are the freezing scan (probably partly liquid). The others are
      mid-run and recover at the next temperature; in them Mustrain and DisplaceX jump
      (DisplaceX up to 114 vs typical ±5), so the fit wandered or the pattern changed (spotty ice?)
- [x] Freezing scan (T_freeze) left out of the CTE fit (`CTE_SETTINGS['exclude_freeze_scan']`)
- [x] AFP16 (4 frozen scans, one bad) and WT20 (2) excluded from the sample list in `main.py`
- [ ] Last-cycle shift/esd is still ~0.3-1 esd (median 0.7) when GSAS-II stops (chi2 change < 0.1%).
      16:27 run: the parameter is Mustrain (216 of 240 rows) or Scale (23), not a/c, so lattice
      results are unaffected; matters only if microstrain is going to be interpreted.
      Flag threshold set to 1.0 for now. `Last_shift_param` records which parameter moved most;
      if it's always the same one, try a tighter convergence (`min dM/M` 1e-4) or more cycles

## 1. Refinement problems that change results

- [x] **1a** `initial_refine`: duplicate dict keys dropped the axial Mustrain setting (step 2) and the
      DisplaceX/Background step (step 5). Rewritten with one key per step.
- [x] **1b** Stop refining Zero on the sample. In 6.5-18.5° 2θ, Zero and DisplaceX shift peaks almost
      identically (cos2θ only goes 0.994 to 0.948), so refining both makes the split arbitrary.
      Decision: Zero fixed at 0, DisplaceX is the only peak-offset parameter.
- [x] **1c** DisplaceX refined per temperature together with Dij in `additional_seq_refine`.
      Tested 2026-10-09 (runs 15:53 refined vs 16:03 fixed): per-temperature DisplaceX scatters
      less than its esd, shows no drift, and is anti-correlated with the a/c scatter. Fixing it
      halved the scatter of a and c around their trends and the CTE uncertainty (median ratio
      0.43 for a, 0.53 for c), and lowered CTEs by ~1 sigma, more for AFP than WT (the refined
      run made AFP look ~7e-6/K higher than WT in a; fixed: ~2e-6/K).
      Decision: `REFINE_DISPLACEMENT_PER_T = False` (`config.py`); DisplaceX is refined only in
      `initial_refine`. Note DisplaceX also shifts absolute a/c (median DisplaceX 42 µm fixed vs
      5 µm refined moved a by 0.006 Å), consistent with the distance-scale caveat below
- [x] **1d** `calculate_cte`: weighted fit (`np.polyfit(..., w=1/sigma, cov=True)`), optional
      `exclude_flagged=True` to leave out `Fit_OK == False` points, skips unrefined temperatures.
      Also reports `Red_chi2_fit` (scatter around the line vs error bars)
- [ ] **1e** (later) Glass background: subtraction is commented out. The glass capillary sits in the
      beam and adds most of the background. Subtracting the glass image before integration has
      made GSAS-II refinements fail in the past, so for now the background function has to carry
      it; check that 5 terms are enough (low Durbin-Watson = not enough)

## Structure model

- [~] `ice.cif` has U_iso = 0 for every atom (no thermal motion), so calculated intensities fall
      off too slowly with angle. Implemented: `UISO_SETTINGS` (`config.py`). Default `'refine'`
      sets O 0.02 / H 0.04 Å² as a start, refines O's U_iso with the cell on the first scan,
      then holds it; `Uiso_O` column in the results. To do: run with `'refine'` and with
      `'cif'` (old behavior) and compare Rwp, a, c and CTE
- [ ] Peak width / crystallite size (goal: compare FWHM or size between AFP and WT).
      Refining Size makes the fit unstable because the data can't separate it from other broadening:
  - Size broadening goes as 1/cos(theta), which changes by only ~1% over 6.5-18.5 deg 2theta, so it
    is indistinguishable from constant-width terms (instrument W/X, and the geometric broadening
    of a ~mm larva at 150 mm). Microstrain (tan theta, varies 3x over the range) is separable
  - Instrument FWHM (x.instprm) is ~0.18 deg at 12 deg 2theta. Scherrer broadening is ~0.03 deg
    for 100 nm crystallites and ~0.003 deg for 1 um, so only crystallites below ~100-200 nm
    would show up at all. Ice grains in frozen tissue are typically microns (rings are spotty)
  - Options, in order of usefulness:
    1. [~] Ring spottiness from the 2D images (intensity variation along each ring vs azimuth):
       measures the number of diffracting grains, i.e. grain size in the micron range where
       AFP's recrystallization inhibition would act. Implemented: `python -m scripts.spottiness`
       (settings in `SPOTTINESS_SETTINGS`); stats and report pick it up automatically.
       To do: run on the real images, check the numbers look sensible (smooth vs spotty
       images by eye), then read the AFP vs WT test
    2. Empirical peak width: FWHM of a few strong reflections per scan (from the refined
       profile or single-peak fits), compared between groups as excess over the LaB6 width.
       Larva size/position adds geometric width, so compare with care
    3. Keep Size fixed (large) and use Mustrain as the width parameter, as now
- [ ] Ag Kalpha1/Kalpha2 doublet: `x.instprm` uses one averaged wavelength (0.56083). The doublet
      splitting grows from 0.05 to 0.15 deg 2theta across the range, which is folded into U,V,W
      from the LaB6 fit. Fine for positions; for width analysis, refit the LaB6 instrument file
      with Lam1/Lam2 (in GSAS-II, saved as a new file outside `controls/`)
      Plan (Mauricio): recalibrate the LaB6 standard by hand in GSAS-II, fitting the Ka2/Ka1
      ratio. Matters mostly at the high-angle end of this range
- [ ] Preferred orientation: test spherical harmonics (order 2-4) on the first scan; spotty
      tissue ice may carry texture that the near-full-ring integration doesn't average out
- [x] Refining O and H positions: not worthwhile with this data (d >= 1.75 Å, ~a dozen
      reflections, H nearly invisible to X-rays, intensities affected by spotty ice). Positions
      don't change a and c. Keep the CIF positions fixed

## Calibration / sample geometry

- [ ] Sample-to-detector distance: larva at the capillary tip is not where the LaB6 was. This scales
      a and c together and can't be fixed by Zero or DisplaceX. CTE is unaffected if the distance is
      constant; absolute a, c are not.
  - [ ] Compare a, c at one temperature with literature ice Ih (Röttger et al., 1994)
  - [ ] If possible, measure a standard in the same capillary-tip position for future runs
- [x] Polarization is 0.99 in `.imctrl` and `x.instprm`. Decision: keep it, per the lab source
      technician (it affects intensities/Rwp, not peak positions)
- [ ] Larvae can move or dehydrate during a run, and ice in tissue can be spotty (few large crystals).
      Watch for single-temperature Rwp jumps and check those images by eye

## Statistics

- [x] Record per sample: temperature range measured (start, end, step) and the freezing point
      bracket (`temperature_info` → *Samples* sheet and `output/sample_info.csv`). Also counts
      scans rejected as unfrozen below the freezing point, which should be 0
- [x] Statistical tests (`python -m scripts.stats`, 2026-10-09 16:27 run). One value per larva;
      permutation tests on the difference in means (AFP vs WT shuffled within start temperature),
      bootstrap 95% CI, Mann-Whitney cross-check, Holm correction over the 9 tests.
  - Collection: each sample starts at 260 K or 280 K and is cooled in 1 K steps, scanned at each
    step. So this is a freezing (supercooling) temperature, known to within 1 K
  - Start temperature: 280 K starts froze 2.4 K warmer (95% CI 1.0-3.7 K, p = 0.005,
    Holm 0.04; AFP +2.2 K, WT +2.6 K separately, p ~0.05 each)
  - Caveat: start temperature is confounded with collection batch (samples 2-8 started at 260 K,
    15-24 at 280 K), so a batch/cohort difference would look the same
  - AFP vs WT freezing: AFP 1.5 K colder (CI 0.3-2.7 K, p = 0.06, Holm 0.38): suggestive only
  - AFP vs WT CTE, a, c, c/a at 248 K: no difference. CTE a +1.1 (CI -1.3 to 3.4) e-6/K,
    CTE c -0.1 (CI -4.0 to 3.5) e-6/K, c/a -4e-5 (CI -1.5e-4 to 6e-5)
- [x] Present: HTML report (`python -m scripts.report` → `output/report.html`; runs the stats first)

## 2. Code that will crash or misbehave

- [x] `analysis.py` reads `seq_results.xlsx`, which `save_data` never writes (point it at the csv)
- [x] `analysis.py` runs plots at import; wrap in `if __name__ == "__main__":`
- [x] `config.py`: create `projects/` and `output/` if missing (`mkdir(exist_ok=True)`)
- [x] `config.py`: give project files a `.gpx` extension
- [x] `main.py`: `try/except` around each sample so one failure doesn't lose the whole run
- [x] `integrate_images`: check the list returned by `Integrate()` before taking `[0]`
- [x] `requirements.txt`: add `pandas` and `openpyxl`; re-save as UTF-8 (added unpinned: re-run
      `pip freeze` on the local machine to pin the installed versions)

## 3. Things that could quietly be wrong

- [ ] `remove_unfrozen`: fixed intensity cutoff depends on exposure; consider a peak-to-median test.
      Check which temperatures it drops per sample
- [ ] `remove_unfrozen`: use `get_sample_name()` so PWDR and IMG names match
      (then `remove_orphan_images` is no longer needed)
- [ ] `LoadProfile` replaces the wavelength with `x.instprm`; confirm it matches `.imctrl` (both 0.56083 now)
- [ ] `copy_displaceX`: copy the list (`list(displace_x)`) instead of sharing one object
- [ ] Size limits in `first_seq_refine` do nothing while Size isn't refined

## Process / project setup

- [x] Fix GitHub access (Claude GitHub App on `mceballos44/ScriptingGSAS`) and push the local commits
- [x] `CLAUDE.md` with working rules for Claude sessions
- [x] `README.md` documenting what each part does
- [ ] Move parameters out of code into one config file: sample list (`main.py`), paths and
      file names (`config.py`), refinement settings and thresholds. Sample metadata is the same
      for every sample, so no per-sample metadata files; the only per-sample difference is the
      temperature range, which can be read from the image file names (record it in the results)
- [ ] Enforce read-only `data/` and `controls/` with deny rules in `.claude/settings.json`
- [ ] On the local machine: run `/fewer-permission-prompts` after a few sessions
- [x] Quieter terminal: GSAS-II output goes to `output/logs/<sample>.log`, terminal shows a
      summary line per sample and GSAS-II warnings/errors (`QUIET_GSAS` in `config.py`)
- [ ] HTML report per run: Rwp and a/c vs T with flags, CTE table, worst fits
- [ ] Once the pipeline is stable: package "run a new sample" as a skill (`skill-creator`)

## 4. Cleanup

- [ ] Remove unused `columns` list in `results.py`
- [ ] Remove `background_list = 0` in `config.py`
- [ ] `ice_cif` defined in four files; keep one
- [ ] Empty `stats_summary` in `analysis.py`
- [ ] Duplicate imports and outdated header comment in `main.py`
