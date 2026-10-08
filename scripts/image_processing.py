# Importing images and integrating them into 1d patterns
# 
# Includes loading set of images, copying controls and mask, integrating entire set

### Written by Mauricio Ceballos, Joester Group, Northwestern University
### Last updated: 2026-09-22

import re
import glob
from zipfile import ZipFile
import G2script as G2sc
from pathlib import Path
import numpy as np
from scripts.config import setup,BASE_DIR,DATA_DIR,PROJECT_DIR,OUTPUT_DIR,CONTROLS_DIR
from scripts.config import BACKGROUND_FILE,CONTROLS_FILE,MASK_FILE,INSTRUMENT_FILE

def integrate_images(
    gpx=None,
    samples=None
):
    """
    Loads a single image into GSAS,
    sets image controls, mask, and background.
    Integrates the image and assigns histogram to phase,
    saves the project.

    Prepares for next step: refinement of the integrated data.
    """
    # Check project
    if gpx is None:
        raise RuntimeError("No GSAS project was specified, check gpx=...")
    # Samples directory
    sample_loc = DATA_DIR / samples

    def get_temperature(filepath):
        match = re.search(r'_(\d+)$',filepath.stem)
        if match is None:
            raise ValueError(
                f"Could not fine temperature in filename: {filepath.name}"
            )
        temperature = float(match.group(1))
        return temperature

    # Find cbf files
    image_files = sorted(sample_loc.glob('*.cbf'))
    # Sort numerically
    image_files = sorted(image_files, key=get_temperature)
    
    
    for image_file in image_files:
        temperature = get_temperature(image_file)
        # Troubleshooting prints
        #print(f"Loading {image_file.name}")
        #print(f"Temperature: {temperature} K")

        image_list = gpx.add_image(
            str(image_file),
            fmthint="CBF"
        )

        if not image_list:
            raise RuntimeError(
                f"Failed to import image: {image_file}"
            )
        
        image = image_list[0]
        # Set image controls and mask
        image.loadControls(str(CONTROLS_FILE))
        image.loadMasks(str(MASK_FILE))

        # Remove sag and change det2theta metadata to 0
        controls = image.data["Image Controls"]
        controls["sag"] = 0.0
        controls["det2theta"] = 0

        # Set background sub
        # image.setControlFile(
        #     "background image",
        #     background,
        #     mult=-1.0
        # )

        # Integrate image
        pwdr = image.Integrate()[0]

        if not pwdr:
            raise RuntimeError(
                f"Integration failed for: {image_file.name}"
            )

        # Set temperature
        pwdr.SampleParameters['Temperature'] = temperature
        # Load Instrument parameters
        pwdr.LoadProfile(str(INSTRUMENT_FILE))
        # print(pwdr.name)
        # print(pwdr.SampleParameters)


    gpx.save()
    #print(f"\n\nProject file name: {gpx}")
    return gpx

def remove_unfrozen(gpx,threshold=0.80):
    """
    This function will remove unfrozen image/pwdr data from project
    Maybe check histogram data to check intensities, removes those that do not fit
    """
    unfrozen = []
    # Find unfrozen
    for hist in gpx.histograms("PWDR "):
        # Pwdr data
        y = hist.data["data"][1][1]
        high_intensity = np.percentile(y,99.9)
        if high_intensity < threshold:
            sample_name = hist.name.replace("PWDR ","")
            sample_name = Path(sample_name).stem

            unfrozen.append(sample_name)

            print(
                f"Marked unfrozen: {sample_name}"
            )
    
    # Find matching PWDR + IMG
    keys_to_delete = []

    for key in list(gpx.data.keys()):
        if not (
            key.startswith("PWDR ")
            or key.startswith("IMG ")
        ):
            continue

        # Remove GSAS prefix
        sample_name = key.replace("PWDR ", "")
        sample_name = sample_name.replace("IMG ","")

        # Remove .cbf
        sample_name = Path(sample_name).stem

        if sample_name in unfrozen:
            keys_to_delete.append(key)
    # Delete
    for key in keys_to_delete:
        print(f"Deleting: {key}")
        # Delete data
        del gpx.data[key]
        # Delete tree
        gpx.names = [
            item for item in gpx.names
            if item[0] != key
        ]
    # Refresh GSAS IDs
    gpx.update_ids()
    gpx.save()
    return unfrozen

def get_sample_name(name):

    name = name.replace("PWDR ", "")
    name = name.replace("IMG ", "")

    # Remove GSAS integration suffix
    name = name.split(" Azm=")[0]

    # Remove .cbf
    name = Path(name).stem

    return name

def remove_orphan_images(gpx):

    # Get sample names that still have PWDR data
    pwdr_samples = set()

    for hist in gpx.histograms("PWDR "):
        sample_name = get_sample_name(hist.name)
        pwdr_samples.add(sample_name)

    # Find images without corresponding PWDR data
    images_to_delete = []

    for key in list(gpx.data.keys()):

        if not key.startswith("IMG "):
            continue

        sample_name = get_sample_name(key)

        if sample_name not in pwdr_samples:
            images_to_delete.append(key)

    # Delete orphan images
    for key in images_to_delete:

        print(f"Deleting orphan image: {key}")

        del gpx.data[key]

        gpx.names = [
            item for item in gpx.names
            if item[0] != key
        ]

    gpx.update_ids()
    gpx.save()
    return gpx

def assign_phase_all(gpx, phase_name):
    """
    Assigns existing phase to all pwdr data
    Inputs:
    gpx: gpx from setup()
    phase_name: (str)
    """
    phase = gpx.phase(phase_name)

    # Get only PWDR histograms
    pwdr_histograms = [
        hist for hist in gpx.histograms()
        if hist.name.startswith("PWDR ")
    ]

    # Skip the first PWDR histogram, first contains refined data
    for hist in pwdr_histograms[1:]:

        gpx.link_histogram_phase(hist, phase)

        print(f"Linked {phase_name} -> {hist.name}")

    gpx.save()
    return gpx

def assign_phase_one(gpx, phase_name):
    """
    Assigns existing phase to all pwdr data
    Inputs:
    gpx: gpx from setup()
    phase_name: (str)
    """
    phase = gpx.phase(phase_name)

    hist = gpx.histogram(0)

    gpx.link_histogram_phase(hist, phase)

    print(f"Linked {phase_name} -> {hist.name}")
    gpx.save()
    return gpx


#Testing the single image function
# from scripts.config import setup
#gpx = setup()
# print(f'BKG DIR: \n{BASE_DIR}\n')
##integrate_images(gpx=gpx,samples ='AFP/AFP3')
# remove_unfrozen(gpx=gpx)
# remove_orphan_images(gpx=gpx)
# assign_phase_one(gpx=gpx,phase_name='ice')