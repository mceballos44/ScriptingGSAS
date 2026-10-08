from pathlib import Path
import G2script as G2sc
import scripts.config
# Overall flow:
# We start with setup.py, get the settings correct there
# 
# setup.py sets directories, makes initial project with
# background image and ice phase
# 
# image_processing contains a couple functions
# the first adds an image to the initial project, sets bkgrnd img,
# mask, control, integrates
#
# refinement.py contains a function that performs basic refinement
# of the first image
# -----------------------------
# Project paths
# -----------------------------
import os
import G2script as G2sc
from pathlib import Path
import pandas as pd
# Import directories
from scripts.config import setup,BASE_DIR,DATA_DIR,PROJECT_DIR,OUTPUT_DIR,CONTROLS_DIR
from scripts.config import BACKGROUND_FILE,CONTROLS_FILE,MASK_FILE,INSTRUMENT_FILE
from scripts.results import extract_data, combine_data, save_data
import scripts.image_processing as ip
import scripts.refinement as rf
# Will always use this cif file for the phase
ice_cif = DATA_DIR / "ice.cif"

def full_analysis(sample_name):
    """
    Run the complete GSASII analysis for a single sample
    """
    
    gpx = setup(sample_name=sample_name)
    
    # Determine sample type:
    if sample_name.startswith('AFP'):
        sample_folder = f'AFP/{sample_name}'
    elif sample_name.startswith('WT'):
        sample_folder = f'WT/{sample_name}'
    else:
        raise ValueError(
            f'Unknown sample type: {sample_name}'
        )
    ip.integrate_images(gpx=gpx,samples=sample_folder)
    ip.remove_unfrozen(gpx=gpx)
    ip.remove_orphan_images(gpx=gpx)

    ip.assign_phase_one(gpx=gpx,phase_name='ice')

    rf.initial_refine(gpx=gpx)
    rf.first_seq_refine(gpx=gpx)
    rf.additional_seq_refine(gpx=gpx)
    df = extract_data(gpx=gpx,sample_name=sample_name)
    
    return df
# -----------------------------
# Main workflow
# -----------------------------

samples = [
    "AFP3",
    "AFP4",
    "AFP5",
    "AFP6",
    "AFP7",
    "AFP8",
    "AFP15",
    "AFP16",
    "AFP17",
    "AFP18",
    "AFP21",
    "WT2",
    "WT3",
    "WT4",
    "WT5",
    "WT6",
    "WT7",
    "WT20",
    "WT21",
    "WT22",
    "WT23",
    "WT24",
]

all_results = []

for sample_name in samples:

    print()
    print(f"Running {sample_name}")
    print()

    df = full_analysis(sample_name)

    all_results.append(df)


full_df = combine_data(all_results)

save_data(full_df)

print(full_df)