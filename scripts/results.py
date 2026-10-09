# Functions here will extract parameters from sequential refinement
# Also set up nice tables
# We want the data to be setup like this
# Sample, T, A, sigmaA, C, sigmaC, Diso, Dsigma, strain,sigmaS, Rwp
import os
import G2script as G2sc
from pathlib import Path
import numpy as np
import pandas as pd
from datetime import datetime
# Import directories
from scripts.config import setup,BASE_DIR,DATA_DIR,PROJECT_DIR,OUTPUT_DIR,CONTROLS_DIR
from scripts.config import BACKGROUND_FILE,CONTROLS_FILE,MASK_FILE,INSTRUMENT_FILE
# Will always use this cif file for the phase
ice_cif = DATA_DIR / "ice.cif"

columns=['Sample','T','A','sigmaA','C','sigmaC',
                 'D11','sigma_D11','D33','sigma_D33',
                 'Astrain','sigma_Astrain','Cstrain',
                 'sigma_Cstrain','Rwp']

def _value_esd(seq, x, var):
    """
    Safe wrapper around seq.get_Variable: returns (value, esd) with NaN
    when the variable is missing or was not refined
    """
    result = seq.get_Variable(x, var)
    if result is None:
        return np.nan, np.nan
    val, esd = result
    if esd is None:
        esd = np.nan
    return val, esd

def _profile_stats(hist_data, rwp_gsas):
    """
    Rp and Durbin-Watson computed from the observed/calculated pattern
    arrays of one histogram. GSAS-II only stores these residuals for the
    last single-pattern refinement, not for every sequential step, so
    they are recomputed here.

    The Rwp from the same arrays is compared with GSAS-II's own Rwp; if
    they disagree the arrays are stale and NaN is returned instead.
    """
    nan = {'Rp': np.nan, 'Durbin_Watson': np.nan}
    try:
        x, yo, w, yc = (np.ma.asarray(a) for a in hist_data['data'][1][:4])
    except (KeyError, IndexError, TypeError):
        return nan
    lo, hi = hist_data['Limits'][1]
    use = (~np.ma.getmaskarray(x)) & (x >= lo) & (x <= hi) & (w > 0)
    yo, w, yc = (np.asarray(a)[use] for a in (yo, w, yc))
    if len(yo) < 3 or not np.any(yc):
        return nan
    diff = yo - yc
    rwp = 100 * np.sqrt(np.sum(w * diff**2) / np.sum(w * yo**2))
    if not np.isfinite(rwp_gsas) or abs(rwp - rwp_gsas) > 0.05 * rwp_gsas:
        return nan
    wdiff = np.sqrt(w) * diff
    return {
        'Rp': 100 * np.sum(np.abs(diff)) / np.sum(yo),
        'Durbin_Watson': np.sum(np.diff(wdiff)**2) / np.sum(wdiff**2),
    }

def _last_shift_esd(seq_results):
    """
    Largest |shift/esd| in the final least-squares cycle, and the
    parameter it belongs to. GSAS-II's own 'Max shft/sig' is the total
    change from the starting values, which is large whenever a parameter
    moved a lot, even in a converged fit.
    """
    shifts = seq_results.get('Rvals', {}).get('lastShifts')
    if not shifts:
        return np.nan, ''
    sig = dict(zip(seq_results.get('varyList', []), seq_results.get('sig', [])))
    ratios = {
        name: abs(shift / sig[name]) for name, shift in shifts.items()
        if sig.get(name) not in (None, 0) and np.isfinite(sig[name])
    }
    if not ratios:
        return np.nan, ''
    worst = max(ratios, key=ratios.get)
    return ratios[worst], worst

def _fit_metrics(seq_results, hist_data):
    """
    Pull the goodness-of-fit numbers for one histogram.

    seq_results['Rvals'] is written by the least-squares engine
    (GSASIIstrMain.RefineCore); Rp and Durbin-Watson are computed from
    the pattern arrays (see _profile_stats).

    Note: GOF assumes intensities are in counts. These integrated
    intensities are much smaller than counts, so GOF comes out far below
    1; compare it between temperatures, not with 1.
    """
    rvals = seq_results.get('Rvals', {})
    rwp = rvals.get('Rwp', np.nan)
    gof = rvals.get('GOF', np.nan)
    last_shift, last_shift_param = _last_shift_esd(seq_results)
    metrics = {
        'Rwp': rwp,                                     # weighted profile R (%)
        'GOF': gof,                                     # relative use only, see above
        'Red_chi2': gof**2,
        'Nobs': rvals.get('Nobs', np.nan),
        'Nvars': rvals.get('Nvars', np.nan),
        'Converged': rvals.get('converged', None),
        'DelChi2': rvals.get('DelChi2', np.nan),        # last relative change in chi2
        'Last_shift_esd': last_shift,                    # final-cycle shift/esd
        'Last_shift_param': last_shift_param,           # parameter with that shift
        'Total_shift_esd': rvals.get('Max shft/sig', np.nan),  # change from start
        'SVD_singular': rvals.get('SVD0', 0),           # >0 means correlated/undetermined params
        'Aborted': rvals.get('Aborted', False),
        'Refine_msg': rvals.get('msg', '').strip(),
    }
    metrics.update(_profile_stats(hist_data, rwp))      # Rp, Durbin_Watson
    return metrics

def extract_data(gpx,sample_name):
    """
    This function should extract the data from the current gpx and 
    store data in dataframe
    """
    rows = []
    seq = gpx.seqref()

    # print(f'\n\nParamlist: {seq.get_ParmList(0)}')
    # print(f'\n\nCell and ESD format: {seq.get_cell_and_esd('ice',0)}')
    # We can use dictionaries instead of constantly appending data
    for x,hist_name in enumerate(seq.histograms()):
        # Sequential refinement stops at the first failed histogram,
        # anything after that has no results
        if hist_name not in seq.data:
            print(f"No sequential results for {hist_name}, refinement stopped before it")
            rows.append({
                'Sample': sample_name,
                'Seq_index': x,
                'Histogram': hist_name,
                'T': gpx.histogram(hist_name).SampleParameters['Temperature'],
                'Refined': False,
            })
            continue
        cell, cellESD, _ = seq.get_cell_and_esd('ice',x)
        ref_data = seq.RefData(x)
        seq_results = ref_data[0]
        D11, sigma_D11 = _value_esd(seq, x, f'0:{x}:D11')
        D33, sigma_D33 = _value_esd(seq, x, f'0:{x}:D33')
        Mustrain_a, sigma_Mustrain_a = _value_esd(seq, x, f'0:{x}:Mustrain;a')
        Mustrain_i, sigma_Mustrain_i = _value_esd(seq, x, f'0:{x}:Mustrain;i')
        # Sample displacement (mm): tracks the larva moving between temperatures
        DisplaceX, sigma_DisplaceX = _value_esd(seq, x, f':{x}:DisplaceX')
        row = {
            'Sample': sample_name,
            'Seq_index': x,
            'Histogram': hist_name,
            'T': seq_results['parmDict'][f':{x}:Temperature'],
            'Refined': True,
            'A': cell[0],
            'sigmaA': cellESD[0],
            'C': cell[2],
            'sigmaC': cellESD[2],
            'D11': D11,
            'sigma_D11': sigma_D11,
            'D33': D33,
            'sigma_D33': sigma_D33,
            'Mustrain_a': Mustrain_a,
            'sigma_Mustrain_a': sigma_Mustrain_a,
            'Mustrain_i': Mustrain_i,
            'sigma_Mustrain_i': sigma_Mustrain_i,
            'DisplaceX': DisplaceX,
            'sigma_DisplaceX': sigma_DisplaceX,
        }
        row.update(_fit_metrics(seq_results, ref_data[1]))
        rows.append(row)
    df = pd.DataFrame(rows)
    return df

def temperature_info(sample_name, all_temps, frozen_temps):
    """
    Summarize the temperature run of one sample and where it froze.

    Samples are cooled in steps, so the freezing temperature is bracketed
    by the warmest scan with ice (T_freeze) and the scan just before it
    with no ice (T_last_unfrozen).

    Inputs:
    all_temps: temperatures of every integrated image
    frozen_temps: temperatures left after remove_unfrozen
    """
    all_temps = sorted(all_temps)
    frozen_temps = sorted(frozen_temps)
    unfrozen_temps = [t for t in all_temps if t not in frozen_temps]

    t_freeze = max(frozen_temps) if frozen_temps else np.nan
    above = [t for t in unfrozen_temps if t > t_freeze]
    below = [t for t in unfrozen_temps if t < t_freeze]
    steps = np.diff(all_temps)

    return {
        'Sample': sample_name,
        'Group': ''.join(c for c in sample_name if c.isalpha()),
        'T_start': max(all_temps) if all_temps else np.nan,  # cooling run: starts warm
        'T_end': min(all_temps) if all_temps else np.nan,
        'T_step': float(np.median(steps)) if len(steps) else np.nan,
        'N_images': len(all_temps),
        'N_frozen': len(frozen_temps),
        'T_freeze': t_freeze,                               # warmest scan with ice
        'T_last_unfrozen': min(above) if above else np.nan, # scan just before freezing
        # Scans rejected as unfrozen below the freezing point: should be 0,
        # otherwise the intensity cutoff dropped frozen frames (check them)
        'N_unfrozen_below_freeze': len(below),
        'T_unfrozen_below_freeze': ', '.join(f'{t:g}' for t in below),
    }

# ---------------------------------
# Fit quality checks
# ---------------------------------

def _robust_high(values, n_mad):
    """
    True where a value is unusually high compared to the rest of the sample:
    above median + n_mad * (scaled median absolute deviation)
    """
    median = values.median()
    mad = 1.4826 * (values - median).abs().median()
    if not np.isfinite(mad) or mad == 0:
        return pd.Series(False, index=values.index)
    return values > median + n_mad * mad

def _lattice_outlier(sample_data, parameter, n_sigma):
    """
    Distance of each point from a straight-line fit of parameter vs T,
    in units of a robust spread of the residuals. Points that sit off the
    line are suspicious since CTE assumes linear behavior.
    """
    good = sample_data[[parameter, 'T']].dropna()
    z = pd.Series(np.nan, index=sample_data.index)
    if len(good) < 4:
        return z
    slope, intercept = np.polyfit(good['T'], good[parameter], 1)
    resid = good[parameter] - (slope * good['T'] + intercept)
    spread = 1.4826 * (resid - resid.median()).abs().median()
    if spread == 0:
        return z
    z[good.index] = resid / spread
    return z

def flag_fit_quality(
    df,
    rwp_tol=0.0,
    n_mad=3.0,
    max_shift_esd=0.1,
    lattice_sigma=3.0
):
    """
    Add quality flags to the results table. Flags are computed within each
    sample, in sequential refinement order (each histogram starts from the
    previous one's result, so a worse fit than the step before means the
    refinement got harder or went wrong at that temperature).

    Inputs:
    df: dataframe from extract_data / combine_data
    rwp_tol: relative Rwp rise allowed before flagging (0.02 = 2%).
        0 flags any increase.
    n_mad: how many robust deviations above the sample median counts
        as an Rwp or GOF outlier
    max_shift_esd: largest acceptable final-cycle shift/esd (converged fits are < 0.1)
    Durbin-Watson is flagged when unusually low for the sample (n_mad robust
    deviations below its median): residuals more correlated than usual,
    i.e. a systematic misfit of peak shapes/background. It is not compared
    with a fixed value because finely sampled patterns sit well below 2.
    lattice_sigma: robust z-score for A or C to count as off the linear trend

    New columns:
    Rwp_change: Rwp minus Rwp of the previous temperature in the sequence
    Rwp_increase: True where Rwp went up from the previous temperature
    Rwp_outlier, GOF_outlier: much worse than typical for that sample
    DW_outlier: Durbin-Watson much lower than typical for that sample
    A_trend_z, C_trend_z: distance from a linear A(T) / C(T) fit
    Fit_flags: text list of every check that failed
    Fit_OK: True when no check failed
    """
    df = df.copy()
    if 'Refined' not in df:
        df['Refined'] = True
    df = df.sort_values(['Sample', 'Seq_index'])

    for col in ['Rwp_change', 'A_trend_z', 'C_trend_z']:
        df[col] = np.nan
    for col in ['Rwp_increase', 'Rwp_outlier', 'GOF_outlier', 'DW_outlier']:
        df[col] = False

    for sample, sample_data in df.groupby('Sample', sort=False):
        refined = sample_data[sample_data['Refined'] == True]
        idx = refined.index

        rwp_change = refined['Rwp'].diff()
        df.loc[idx, 'Rwp_change'] = rwp_change
        df.loc[idx, 'Rwp_increase'] = (
            refined['Rwp'] > refined['Rwp'].shift() * (1 + rwp_tol)
        )
        df.loc[idx, 'Rwp_outlier'] = _robust_high(refined['Rwp'], n_mad)
        df.loc[idx, 'GOF_outlier'] = _robust_high(refined['GOF'], n_mad)
        if 'Durbin_Watson' in refined:
            # Low DW is bad: flip the sign and reuse the high-outlier test
            df.loc[idx, 'DW_outlier'] = _robust_high(-refined['Durbin_Watson'], n_mad)
        df.loc[idx, 'A_trend_z'] = _lattice_outlier(refined, 'A', lattice_sigma)
        df.loc[idx, 'C_trend_z'] = _lattice_outlier(refined, 'C', lattice_sigma)

    def reasons(row):
        if row['Refined'] != True:
            return 'not refined'
        flags = []
        if row['Rwp_increase']:
            flags.append(f"Rwp up {row['Rwp_change']:+.2f}")
        if row['Rwp_outlier']:
            flags.append('Rwp outlier')
        if row['GOF_outlier']:
            flags.append('GOF outlier')
        if row.get('Converged') is not None and row.get('Converged') == False:
            flags.append('not converged')
        if row.get('Aborted') == True:
            flags.append('aborted')
        if row.get('Last_shift_esd', 0) > max_shift_esd:
            flags.append('large shift/esd')
        if row.get('SVD_singular', 0) > 0:
            flags.append('singular params')
        if row['DW_outlier']:
            flags.append('correlated residuals')
        if abs(row['A_trend_z']) > lattice_sigma:
            flags.append('A off trend')
        if abs(row['C_trend_z']) > lattice_sigma:
            flags.append('C off trend')
        return '; '.join(flags)

    df['Fit_flags'] = df.apply(reasons, axis=1)
    df['Fit_OK'] = df['Fit_flags'] == ''
    return df

def quality_summary(df):
    """
    One row per sample: Rwp range, how many temperatures were flagged
    and which ones
    """
    rows = []
    for sample, sample_data in df.groupby('Sample', sort=False):
        flagged = sample_data[~sample_data['Fit_OK']]
        increased = sample_data[sample_data['Rwp_increase'] == True]
        rows.append({
            'Sample': sample,
            'N_temps': len(sample_data),
            'N_refined': int((sample_data['Refined'] == True).sum()),
            'Rwp_min': sample_data['Rwp'].min(),
            'Rwp_median': sample_data['Rwp'].median(),
            'Rwp_max': sample_data['Rwp'].max(),
            'GOF_median': sample_data['GOF'].median(),
            'N_Rwp_increase': len(increased),
            'T_Rwp_increase': ', '.join(f'{t:g}' for t in increased['T']),
            'N_flagged': len(flagged),
            'T_flagged': ', '.join(f'{t:g}' for t in flagged['T']),
        })
    return pd.DataFrame(rows)

def combine_data(sample_list):
    full_df = pd.concat(
        sample_list,
        ignore_index=True
    )
    
    # Save dataframe as readable csv/excel file
    return full_df

def save_data(df, sample_info=None):
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")
    excel_file = OUTPUT_DIR / f'seq_results_{timestamp}.xlsx'
    csv_file = OUTPUT_DIR / 'seq_results.csv'
    
    with pd.ExcelWriter(excel_file) as writer:
        df.to_excel(writer, sheet_name='Results', index=False)
        if 'Fit_OK' in df:
            quality_summary(df).to_excel(
                writer, sheet_name='Fit_quality', index=False
            )
        if sample_info is not None:
            sample_info.to_excel(writer, sheet_name='Samples', index=False)
    df.to_csv(
        csv_file,
        index=False
    )
    print(f"Saved results to:")
    print(excel_file)
    print(csv_file)
    if sample_info is not None:
        info_file = OUTPUT_DIR / 'sample_info.csv'
        sample_info.to_csv(info_file, index=False)
        print(info_file)
    return excel_file
# import scripts.config as sp
# import image_processing as ip
# import refinement as rf
# gpx = sp.setup()
# ip.integrate_images(gpx=gpx,samples ='AFP/AFP3')
# ip.remove_unfrozen(gpx=gpx)
# ip.remove_orphan_images(gpx=gpx)

# ip.assign_phase_one(gpx=gpx,phase_name='ice')

# rf.initial_refine(gpx=gpx)
# rf.first_seq_refine(gpx=gpx)
# rf.additional_seq_refine(gpx=gpx)
# extract_data(gpx=gpx,sample_name='AFP3')
