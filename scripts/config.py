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
    project_file = PROJECT_DIR / sample_name
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