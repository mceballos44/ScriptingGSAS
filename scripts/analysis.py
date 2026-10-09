import os
import G2script as G2sc
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
# Import directories
from scripts.config import setup,BASE_DIR,DATA_DIR,PROJECT_DIR,OUTPUT_DIR,CONTROLS_DIR
#from scripts.config import BACKGROUND_FILE,CONTROLS_FILE,MASK_FILE,INSTRUMENT_FILE
import scripts.image_processing as ip
# Will always use this cif file for the phase
ice_cif = DATA_DIR / "ice.cif"

RESULTS_FILE = OUTPUT_DIR / "seq_results.csv"
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
    df = pd.read_csv(RESULTS_FILE)
    
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

def calculate_cte(parameter, group, sigma_parameter=None, exclude_flagged=False):
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

    sigma_parameter : str or None
        Column with the uncertainty of parameter, used to weight the fit
        (points with larger error bars count less). Defaults to
        "sigma" + parameter, e.g. "sigmaA". If any point is missing an
        uncertainty, that sample falls back to an unweighted fit.

    exclude_flagged : bool
        If True, leave out temperatures with Fit_OK == False
        (see results.flag_fit_quality).

    Returns
    -------
    pandas DataFrame
        Table containing sample, slope, average lattice
        parameter, CTE and its uncertainty, number of points used,
        whether the fit was weighted, and Red_chi2_fit: scatter of the
        points around the line relative to their error bars (about 1 if
        the line and the error bars agree; much larger means extra
        scatter or a curved trend, and Sigma_cte already includes it).
    """
    df = pd.read_csv(RESULTS_FILE)

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

    
    if sigma_parameter is None:
        sigma_parameter = f"sigma{parameter}"

    results = []
    
    # Go through samples one at a time
    for sample in df['Sample'].unique():
        sample_data = df[df['Sample'] == sample].copy()
        n_total = len(sample_data)
        # Temperatures that were never refined have no value
        sample_data = sample_data.dropna(subset=['T', parameter])
        if exclude_flagged and 'Fit_OK' in sample_data:
            sample_data = sample_data[sample_data['Fit_OK'] == True]
        n_used = len(sample_data)

        T = sample_data['T'].to_numpy(dtype=float)
        y = sample_data[parameter].to_numpy(dtype=float)

        # Weight each point by 1/sigma if every point has a usable sigma
        weights = None
        if sigma_parameter in sample_data:
            sigma = sample_data[sigma_parameter].to_numpy(dtype=float)
            if np.all(np.isfinite(sigma)) and np.all(sigma > 0):
                weights = 1.0 / sigma

        slope = sigma_slope = red_chi2 = np.nan
        if n_used >= 4:
            # Calculate slope dL/dT using a (weighted) straight-line fit.
            # cov=True scales the uncertainty by the actual scatter
            # around the line, so underestimated GSAS esds don't make
            # the CTE look more precise than it is
            (slope, intercept), cov = np.polyfit(T, y, 1, w=weights, cov=True)
            sigma_slope = np.sqrt(cov[0, 0])
            resid = y - (slope * T + intercept)
            w2 = 1.0 if weights is None else weights**2
            red_chi2 = np.sum(w2 * resid**2) / (n_used - 2)
        elif n_used >= 2:
            slope = np.polyfit(T, y, 1, w=weights)[0]
            print(f"{sample}: only {n_used} points, no uncertainty on CTE")
        else:
            print(f"{sample}: fewer than 2 usable points, skipped")
        
        # Calculate average lattice parameter
        avg_param = np.mean(y) if n_used else np.nan
        # Calculate CTE
        cte = slope / avg_param
        sigma_cte = sigma_slope / avg_param
        results.append({
            'Sample': sample,
            'Parameter': parameter,
            'Slope': slope,
            'Average': avg_param,
            'CTE': cte,
            'Sigma_cte': sigma_cte,
            'N_points': n_used,
            'N_excluded': n_total - n_used,
            'Weighted': weights is not None,
            'Red_chi2_fit': red_chi2,
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
    df = pd.read_csv(RESULTS_FILE)
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
if __name__ == "__main__":
    plot_rwp(group='AFP')
    plot_parameter('A','sigmaA', group='AFP')
    plot_parameter('C','sigmaC',group='AFP')
    # Sample displacement vs T: a steady drift means the larva moved (TODO 1c)
    plot_parameter('DisplaceX','sigma_DisplaceX',group='AFP')

    a_cte = calculate_cte('A',group='AFP')
    c_cte = calculate_cte('C',group='AFP')

    print(a_cte)
    print(c_cte)