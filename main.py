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
import traceback
import contextlib
import G2script as G2sc
from pathlib import Path
import pandas as pd
# Import directories
from scripts.config import setup,BASE_DIR,DATA_DIR,PROJECT_DIR,OUTPUT_DIR,CONTROLS_DIR
from scripts.config import BACKGROUND_FILE,CONTROLS_FILE,MASK_FILE,INSTRUMENT_FILE
from scripts.results import extract_data, combine_data, save_data
from scripts.results import flag_fit_quality, quality_summary, temperature_info
from scripts.config import QUIET_GSAS, LOG_DIR, GSAS_LOG_KEYWORDS, FIT_CHECKS
from scripts.config import REFINE_DISPLACEMENT_PER_T
from scripts.output_control import gsas_output_to
import scripts.image_processing as ip
import scripts.refinement as rf
# Will always use this cif file for the phase
ice_cif = DATA_DIR / "ice.cif"

def full_analysis(sample_name, refine_displacement=REFINE_DISPLACEMENT_PER_T):
    """
    Run the complete GSASII analysis for a single sample

    refine_displacement: passed to additional_seq_refine; default from
    config.REFINE_DISPLACEMENT_PER_T (False: DisplaceX fixed at the value
    from initial_refine, see TODO 1c)
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
    all_temps = [
        hist.SampleParameters['Temperature']
        for hist in gpx.histograms("PWDR ")
    ]
    ip.remove_unfrozen(gpx=gpx)
    ip.remove_orphan_images(gpx=gpx)
    frozen_temps = [
        hist.SampleParameters['Temperature']
        for hist in gpx.histograms("PWDR ")
    ]
    info = temperature_info(sample_name, all_temps, frozen_temps)

    ip.assign_phase_one(gpx=gpx,phase_name='ice')

    rf.initial_refine(gpx=gpx)
    rf.first_seq_refine(gpx=gpx)
    rf.additional_seq_refine(
        gpx=gpx,
        refine_displacement=refine_displacement
    )
    df = extract_data(gpx=gpx,sample_name=sample_name)
    
    return df, info
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
    # "AFP16",  excluded: only 4 frozen scans, one of them a bad fit
    "AFP17",
    "AFP18",
    "AFP21",
    "WT2",
    "WT3",
    "WT4",
    "WT5",
    "WT6",
    "WT7",
    # "WT20",   excluded: only 2 frozen scans (run ends at 250 K)
    "WT21",
    "WT22",
    "WT23",
    "WT24",
]

all_results = []
all_info = []
failed = {}

for sample_name in samples:

    print(f"Running {sample_name} ...")

    # GSAS-II printouts go to a log file per sample (see config.QUIET_GSAS)
    log_file = LOG_DIR / f"{sample_name}.log"
    if QUIET_GSAS:
        output = gsas_output_to(log_file, GSAS_LOG_KEYWORDS)
    else:
        output = contextlib.nullcontext()

    # One failed sample shouldn't lose the whole run: report it and move on
    try:
        with output:
            df, info = full_analysis(sample_name)
    except Exception as err:
        traceback.print_exc()
        print(f"\n*** {sample_name} failed: {err}")
        if QUIET_GSAS:
            print(f"    GSAS-II output: {log_file}")
        print()
        failed[sample_name] = str(err)
        continue

    all_results.append(df)
    all_info.append(info)

    # One-line summary per sample
    n_refined = int(df['Refined'].sum())
    print(
        f"  {info['T_start']:g} -> {info['T_end']:g} K, froze at {info['T_freeze']:g} K, "
        f"{n_refined}/{len(df)} temperatures refined, "
        f"Rwp {df.get('Rwp', pd.Series(dtype=float)).min():.2f}-"
        f"{df.get('Rwp', pd.Series(dtype=float)).max():.2f}%"
    )
    if info['N_unfrozen_below_freeze']:
        print(
            f"  Check: frames below freezing rejected as unfrozen at "
            f"{info['T_unfrozen_below_freeze']} K"
        )

if not all_results:
    raise RuntimeError("Every sample failed, nothing to save")

full_df = combine_data(all_results)
# Mark temperatures where Rwp went up or other fit checks failed
full_df = flag_fit_quality(full_df, **FIT_CHECKS)

sample_info = pd.DataFrame(all_info)
save_data(full_df, sample_info=sample_info)

# Full table is in output/seq_results.csv; print the per-sample summaries
print()
print(quality_summary(full_df).to_string(index=False))
print()
print(sample_info.to_string(index=False))

if failed:
    print("\nFailed samples:")
    for sample_name, err in failed.items():
        print(f"  {sample_name}: {err}")
