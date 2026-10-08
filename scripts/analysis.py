import os
import G2script as G2sc
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import linregress
# Import directories
from scripts.config import setup,BASE_DIR,DATA_DIR,PROJECT_DIR,OUTPUT_DIR,CONTROLS_DIR
#from scripts.config import BACKGROUND_FILE,CONTROLS_FILE,MASK_FILE,INSTRUMENT_FILE
import scripts.image_processing as ip
# Will always use this cif file for the phase
ice_cif = DATA_DIR / "ice.cif"

RESULTS_FILE = OUTPUT_DIR / "seq_results.xlsx"
# This workbook will contain the functions for rearranging long table data for auto plotting

# First thing I would like to do is get a summary graph of the a and c lattice parameters for all the samples.
# This has to be arranged a certain way:
# Temps 240-280K, samples align automatically to what is available
# Temps - Sample_name - Sample_name -....

# So essentially extend horizontally

def get_parameter_table(df,parameter):
    parameter_table = df.pivot_table(
        index='T',
        columns='Sample',
        values=parameter
    )
    return parameter_table

def plot_parameter(parameter,sigma_parameter,group):
    # Import written excel file
    df = pd.read_excel(RESULTS_FILE)
    
    plt.figure()
    
    # Filter by sample name
    if group is not None:
        df = df[df['Sample'].str.startswith(group)].copy()    
    

    # Extract text group: AFP, WT, etc.
    df["Sample_Group"] = (
        df["Sample"]
        .str.extract(r"([A-Za-z]+)")
    )

    # Extract numeric sample ID: 3, 4, 15, etc.
    df["Sample_Number"] = (
        df["Sample"]
        .str.extract(r"(\d+)")
        .astype(int)
    )

    # Sort by group, sample number, then temperature
    df = df.sort_values(
        ["Sample_Group", "Sample_Number", "T"]
    )
    
    for sample in df["Sample"].unique():
    # Loop through each sample
        # Select rows belonging to this sample
        sample_data = df[
            df['Sample'] == sample
        ].copy()   

        # Plot with error bars if a sigma column was provided
        if sigma_parameter is not None:
            plt.errorbar(
                sample_data['T'],
                sample_data[parameter],
                yerr=sample_data[sigma_parameter],
                label=sample,
                marker='o'
            )
        else:
            plt.plot(
                sample_data['T'],
                sample_data[parameter],
                label=sample,
                marker='o'
            )
    plt.xlabel('Temperature (K)')
    plt.ylabel(parameter)
    plt.legend()
    plt.tight_layout()
    plt.show()
    return

def calculate_cte(parameter,group):
    """
    Calculate linear coefficient of thermal expansion (CTE)
    for each sample.

    CTE = (1 / average lattice parameter) * (slope dL/dT)

    Parameters
    ----------
    df : pandas DataFrame
        Results dataframe.

    parameter : str
        Lattice parameter to use, such as "A" or "C".

    group : str or None
        Optional sample filter.
        Examples:
            "AFP"
            "WT"
        If None, all samples are used.

    Returns
    -------
    pandas DataFrame
        Table containing sample, slope, average lattice
        parameter, and CTE.
    """
    df = pd.read_excel(RESULTS_FILE)

    if group is not None:
        df = df[
            df["Sample"].str.startswith(group)
        ].copy()

    # Extract text group: AFP, WT, etc.
    df["Sample_Group"] = (
        df["Sample"]
        .str.extract(r"([A-Za-z]+)")
    )

    # Extract numeric sample ID: 3, 4, 15, etc.
    df["Sample_Number"] = (
        df["Sample"]
        .str.extract(r"(\d+)")
        .astype(int)
    )

    # Sort by group, sample number, then temperature
    df = df.sort_values(
        ["Sample_Group", "Sample_Number", "T"]
    )

    
    results = []
    
    # Go through samples one at a time
    for sample in df['Sample'].unique():
        sample_data = df[df['Sample'] == sample].copy()
        # Calculate slope dL/dT using linear regression
        fit = linregress(sample_data['T'], sample_data[parameter])
        slope = fit.slope
        sigma_slope = fit.stderr
        
        # Calculate average lattice parameter
        avg_param = sample_data[parameter].mean()
        # Calculate CTE
        cte = slope / avg_param
        sigma_cte = sigma_slope / avg_param
        results.append({
            'Sample': sample,
            'Parameter': parameter,
            'Slope': slope,
            'Average': avg_param,
            'CTE': cte,
            'Sigma_cte': sigma_cte
        })
    
    results_df = pd.DataFrame(results)
    return results_df

def plot_rwp(group):
    """
    Plot Rwp vs temperature for each sample, circling temperatures
    where Rwp went up from the previous step (Rwp_increase) and
    marking any other failed fit check (Fit_OK False) with an x.
    Reads the csv from save_data, which has the flag columns.
    """
    df = pd.read_csv(OUTPUT_DIR / "seq_results.csv")
    if group is not None:
        df = df[df['Sample'].str.startswith(group)].copy()

    plt.figure()
    for sample in df['Sample'].unique():
        sample_data = df[df['Sample'] == sample].sort_values('T')
        line, = plt.plot(
            sample_data['T'],
            sample_data['Rwp'],
            label=sample,
            marker='o'
        )
        increased = sample_data[sample_data['Rwp_increase'] == True]
        plt.scatter(
            increased['T'],
            increased['Rwp'],
            s=150,
            facecolors='none',
            edgecolors=line.get_color()
        )
        bad = sample_data[sample_data['Fit_OK'] == False]
        plt.scatter(
            bad['T'],
            bad['Rwp'],
            marker='x',
            s=80,
            color='black'
        )
    plt.xlabel('Temperature (K)')
    plt.ylabel('Rwp (%)')
    plt.title('Circled: Rwp increase, x: failed fit check')
    plt.legend()
    plt.tight_layout()
    plt.show()
    return

def stats_summary(parameter,group):
    
    return
plot_rwp(group='AFP')
plot_parameter('A','sigmaA', group='AFP')
plot_parameter('C','sigmaC',group='AFP')

a_cte = calculate_cte('A',group='AFP')
c_cte = calculate_cte('C',group='AFP')

print(a_cte)
print(c_cte)