# Statistical comparisons between samples
#
# One value per sample (per larva) is the unit of comparison: CTE of a and c,
# a, c and c/a at a reference temperature, and freezing temperature. Using
# every temperature as a separate data point would overstate the evidence,
# since the temperatures of one larva are not independent.
#
# Tests are permutation tests on the difference in means (exact enough with
# ~10 samples per group and handles the 1 K ties in freezing temperature),
# with a bootstrap 95% interval for the difference and a Mann-Whitney U test
# as a rank-based cross-check. AFP vs WT tests shuffle group labels only
# within the same start temperature (260 or 280 K), so a start-temperature
# effect can't masquerade as a group effect. p-values are also given with a
# Holm correction for the number of tests run together.
#
# Run from the repo root after main.py:  python -m scripts.stats

### Written by Mauricio Ceballos, Joester Group, Northwestern University

import numpy as np
import pandas as pd
from scripts.config import OUTPUT_DIR, CTE_SETTINGS, STATS_SETTINGS
from scripts.analysis import calculate_cte, SAMPLE_INFO_FILE
from scripts.spottiness import SPOTTINESS_SUMMARY_FILE

try:
    from scipy.stats import mannwhitneyu
except ImportError:          # scipy is in requirements.txt, but keep going without it
    mannwhitneyu = None

PER_SAMPLE_FILE = OUTPUT_DIR / "stats_per_sample.csv"
TESTS_FILE = OUTPUT_DIR / "stats_tests.csv"


def per_sample_table(t_ref):
    """
    One row per sample: CTE of a and c (1e-6/K), a, c and c/a at t_ref,
    freezing temperature and start temperature.
    """
    rows = {}
    for parameter in ['A', 'C']:
        fit = calculate_cte(parameter, group=None, t_ref=t_ref, **CTE_SETTINGS)
        for r in fit.itertuples():
            row = rows.setdefault(r.Sample, {'Sample': r.Sample})
            row[f'CTE_{parameter}'] = r.CTE * 1e6
            row[f'sigma_CTE_{parameter}'] = r.Sigma_cte * 1e6
            row[f'{parameter}_Tref'] = r.At_Tref
            row[f'sigma_{parameter}_Tref'] = r.Sigma_at_Tref
            row[f'N_points'] = r.N_points
            row['T_min'], row['T_max'] = r.T_min, r.T_max
    df = pd.DataFrame(rows.values())

    # c/a cancels any error that scales a and c together, such as the
    # sample-to-detector distance, so it is the most robust structural number
    df['C_over_A_Tref'] = df['C_Tref'] / df['A_Tref']
    df['sigma_C_over_A_Tref'] = df['C_over_A_Tref'] * np.sqrt(
        (df['sigma_A_Tref'] / df['A_Tref'])**2
        + (df['sigma_C_Tref'] / df['C_Tref'])**2
    )
    df['Tref_extrapolated'] = (t_ref < df['T_min']) | (t_ref > df['T_max'])

    info = pd.read_csv(SAMPLE_INFO_FILE)
    df = df.merge(info[['Sample', 'Group', 'T_start', 'T_freeze']], on='Sample')

    # Ring spottiness, if scripts/spottiness.py has been run
    if SPOTTINESS_SUMMARY_FILE.exists():
        spots = pd.read_csv(SPOTTINESS_SUMMARY_FILE)
        df = df.merge(spots[['Sample', 'Spottiness', 'Spot_fraction']],
                      on='Sample', how='left')
    return df


def _shuffle_within(labels, strata, rng):
    """Shuffle group labels, but only among samples in the same stratum"""
    shuffled = labels.copy()
    for s in np.unique(strata):
        idx = np.where(strata == s)[0]
        shuffled[idx] = rng.permutation(labels[idx])
    return shuffled


def compare_groups(values, labels, group_a, group_b, strata, n_resamples, rng):
    """
    Compare one quantity between two groups of samples.

    values, labels, strata: arrays with one entry per sample. strata=None
    shuffles labels freely; otherwise only within each stratum.
    Returns means, medians, difference (a - b) with bootstrap 95% interval,
    and permutation and Mann-Whitney p-values.
    """
    keep = np.isfinite(values) & np.isin(labels, [group_a, group_b])
    values, labels = values[keep], labels[keep]
    strata = np.zeros(len(values)) if strata is None else strata[keep]
    xa, xb = values[labels == group_a], values[labels == group_b]

    observed = xa.mean() - xb.mean()

    # Permutation test: how often does a random relabelling give a
    # difference at least as large as the one observed?
    count = 0
    for _ in range(n_resamples):
        shuffled = _shuffle_within(labels, strata, rng)
        diff = values[shuffled == group_a].mean() - values[shuffled == group_b].mean()
        if abs(diff) >= abs(observed) - 1e-12:
            count += 1
    p_perm = (count + 1) / (n_resamples + 1)

    # Bootstrap 95% interval for the difference in means (resample each
    # group within each stratum)
    diffs = np.empty(n_resamples)
    for i in range(n_resamples):
        boot_a, boot_b = [], []
        for s in np.unique(strata):
            in_s = strata == s
            for boot, grp in ((boot_a, group_a), (boot_b, group_b)):
                pool = values[in_s & (labels == grp)]
                if len(pool):
                    boot.append(rng.choice(pool, size=len(pool)))
        diffs[i] = np.concatenate(boot_a).mean() - np.concatenate(boot_b).mean()
    ci_low, ci_high = np.percentile(diffs, [2.5, 97.5])

    p_mwu = np.nan
    if mannwhitneyu is not None and len(xa) and len(xb):
        p_mwu = mannwhitneyu(xa, xb, alternative='two-sided').pvalue

    return {
        'Group_A': group_a, 'Group_B': group_b,
        'N_A': len(xa), 'N_B': len(xb),
        'Mean_A': xa.mean(), 'Mean_B': xb.mean(),
        'Median_A': np.median(xa), 'Median_B': np.median(xb),
        'Diff_A_minus_B': observed, 'CI95_low': ci_low, 'CI95_high': ci_high,
        'p_permutation': p_perm, 'p_mann_whitney': p_mwu,
    }


def holm(pvalues):
    """Holm-Bonferroni adjusted p-values (controls false positives across tests)"""
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (len(p) - rank) * p[i])
        adjusted[i] = min(1.0, running)
    return adjusted


def run_stats(t_ref, n_resamples, seed):
    rng = np.random.default_rng(seed)
    df = per_sample_table(t_ref)
    labels = df['Group'].to_numpy()
    start = df['T_start'].to_numpy()

    tests = []
    # Freezing temperature vs start temperature, within each group
    for group in ['AFP', 'WT']:
        sub = df[df['Group'] == group]
        res = compare_groups(
            sub['T_freeze'].to_numpy(float), sub['T_start'].to_numpy(),
            280, 260, None, n_resamples, rng
        )
        tests.append({'Quantity': 'T_freeze (K)',
                      'Comparison': f'{group}: start 280 vs 260 K', **res})

    # Both groups together: start temperature effect, shuffling only
    # within each group so an AFP/WT difference can't masquerade as it
    res = compare_groups(
        df['T_freeze'].to_numpy(float), start, 280, 260, labels,
        n_resamples, rng
    )
    tests.append({'Quantity': 'T_freeze (K)',
                  'Comparison': 'Both groups: start 280 vs 260 K', **res})

    # AFP vs WT, shuffling only within the same start temperature
    quantities = [
        ('T_freeze', 'T_freeze (K)'),
        ('CTE_A', 'CTE a (1e-6/K)'),
        ('CTE_C', 'CTE c (1e-6/K)'),
        ('A_Tref', f'a at {t_ref:g} K (A)'),
        ('C_Tref', f'c at {t_ref:g} K (A)'),
        ('C_over_A_Tref', f'c/a at {t_ref:g} K'),
    ]
    if 'Spottiness' in df:
        quantities.append(('Spottiness', 'Ring spottiness (CV excess)'))
    for column, name in quantities:
        res = compare_groups(
            df[column].to_numpy(float), labels, 'AFP', 'WT', start,
            n_resamples, rng
        )
        tests.append({'Quantity': name, 'Comparison': 'AFP vs WT', **res})

    tests = pd.DataFrame(tests)
    tests['p_holm'] = holm(tests['p_permutation'])
    return df, tests


if __name__ == "__main__":
    per_sample, tests = run_stats(**STATS_SETTINGS)
    per_sample.to_csv(PER_SAMPLE_FILE, index=False)
    tests.to_csv(TESTS_FILE, index=False)

    pd.set_option('display.width', 200)
    cols = ['Sample', 'Group', 'T_start', 'T_freeze', 'CTE_A', 'sigma_CTE_A',
            'CTE_C', 'sigma_CTE_C', 'A_Tref', 'C_Tref', 'C_over_A_Tref',
            'Tref_extrapolated']
    print(per_sample[cols].round(5).to_string(index=False))
    print()
    show = ['Quantity', 'Comparison', 'N_A', 'N_B', 'Mean_A', 'Mean_B',
            'Diff_A_minus_B', 'CI95_low', 'CI95_high', 'p_permutation',
            'p_mann_whitney', 'p_holm']
    print(tests[show].to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    print(f"\nSaved {PER_SAMPLE_FILE} and {TESTS_FILE}")
