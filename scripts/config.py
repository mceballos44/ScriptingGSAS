from pathlib import Path
import G2script as G2sc

# ---------------------------------
# Setup directories
# ---------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
PROJECT_DIR = BASE_DIR / "projects"
OUTPUT_DIR = BASE_DIR / "output"
CONTROLS_DIR = BASE_DIR / "controls"

BACKGROUND_FILE = DATA_DIR / "glass.cbf"
CONTROLS_FILE = CONTROLS_DIR / ".imctrl"
MASK_FILE = CONTROLS_DIR / ".immask"
INSTRUMENT_FILE = CONTROLS_DIR / "x.instprm"

# Output folders are git-ignored, so create them if this is a fresh checkout
PROJECT_DIR.mkdir(exist_ok=True)
OUTPUT_DIR.mkdir(exist_ok=True)

# ---------------------------------
# Terminal output
# ---------------------------------
# True: GSAS-II's printouts for each sample go to output/logs/<sample>.log
# and the terminal only shows the pipeline summary plus GSAS-II lines
# containing one of the keywords below. False: print everything.
QUIET_GSAS = True
LOG_DIR = OUTPUT_DIR / "logs"
GSAS_LOG_KEYWORDS = ['error', 'warn', 'fail', 'singular', 'abort']

# ---------------------------------
# Refinement settings
# ---------------------------------
# Refine DisplaceX at every temperature in the second sequential pass.
# False (default): keep the value from initial_refine. Per-temperature
# DisplaceX is not determined by the data and trades off with a and c;
# fixing it halved the scatter of a, c and the CTE uncertainty (TODO 1c)
REFINE_DISPLACEMENT_PER_T = False

# ---------------------------------
# CTE fit (analysis.calculate_cte)
# ---------------------------------
CTE_SETTINGS = {
    # Leave out the freezing scan (T_freeze): water is still freezing, so
    # the pattern can be partly liquid / spotty (e.g. AFP5 and AFP8 at 250 K)
    'exclude_freeze_scan': True,
    # Leave out temperatures that failed any fit check (Fit_OK == False)
    'exclude_flagged': False,
    # Individual scans left out of the analysis (the sample stays in):
    # {sample: [T, ...]}. These had Rwp 17-25% in the 2026-10-09 runs, the
    # same with DisplaceX fixed or refined, so the pattern itself is off
    'excluded_scans': {
        'AFP6': [243],
        'AFP17': [246],
        'WT23': [246],
        'WT3': [241, 244, 245],
    },
}

# ---------------------------------
# Statistics (scripts/stats.py)
# ---------------------------------
STATS_SETTINGS = {
    # Temperature (K) at which a, c and c/a are compared between groups.
    # 248 K is inside the measured frozen range of every sample except
    # WT21 (250-255 K, extrapolated 2 K)
    't_ref': 248.0,
    # Random permutations / bootstrap resamples
    'n_resamples': 20000,
    # Seed so the p-values and intervals are the same on every run
    'seed': 0,
}

# ---------------------------------
# Fit quality checks (results.flag_fit_quality)
# ---------------------------------
FIT_CHECKS = {
    # Relative Rwp rise from the previous temperature that gets flagged.
    # First run: rises of 0-3% are noise, real problems were +80% or more
    'rwp_tol': 0.05,
    # Robust deviations from the sample median that count as an outlier
    'n_mad': 3.0,
    # Largest acceptable final-cycle shift/esd. GSAS-II stops when chi2
    # changes by <0.1% per cycle; at that point the last shifts are still
    # ~0.3-1 esd (median 0.7 in the 2026-10-09 runs), so 0.1 flags every
    # row. 1.0 flags only fits still moving by more than their esd
    'max_shift_esd': 1.0,
    # Robust z-score for a or c to count as off the straight-line trend
    'lattice_sigma': 3.0,
}
# ---------------------------------
# Settings, change the quoted pieces
# ---------------------------------
project_name = 'test.gpx'

phase_file = DATA_DIR / 'ice.cif'
background_file = DATA_DIR / 'glass.cbf'
# ---------------------------------
# Setting up base project file
# ---------------------------------
def setup(sample_name):
    project_file = PROJECT_DIR / f"{sample_name}.gpx"
    if project_file.exists():
            project_file.unlink()  # Remove existing project file to avoid conflicts

    gpx = G2sc.G2Project(
                newgpx=project_file,
                author="Mauricio Ceballos, Joester Group, Northwestern University"
            )

    # Add phase to project
    ice_phase = gpx.add_phase(
            str(phase_file),
            phasename="ice",
            fmthint='CIF',
    )

    # Add background file
    background_list = 0
#     background_list = gpx.add_image(
#         str(background_file),
#         fmthint="CBF"
#     )

    gpx.save()
    return gpx