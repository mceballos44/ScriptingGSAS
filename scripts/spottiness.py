# Ring spottiness: how grainy the ice diffraction rings are on the 2D images.
#
# Crystallites of microns or more don't broaden the peaks (see TODO.md, peak
# width), but they do make the rings spotty: with few grains in the beam,
# the intensity along a ring jumps from bin to bin. The size of those jumps,
# after removing what counting noise alone would give, measures how many
# grains contribute, i.e. grain size in the micron range.
#
# For each frozen scan and each reflection in SPOTTINESS_SETTINGS:
#   1. split the ring into azimuth bins (within the integration range)
#   2. net intensity per bin = mean of the ring pixels - mean of nearby
#      background pixels (both from the raw counts)
#   3. CV_obs   = spread of the bin intensities / their mean
#      CV_noise = the spread counting noise alone would give
#      CV_excess = sqrt(CV_obs^2 - CV_noise^2): the grain-statistics part
#   4. N_eff = 1 / CV_excess^2: effective number of grains per bin.
#      Fewer grains (bigger crystals) -> larger CV_excess, smaller N_eff
#   5. Spot fraction: share of bins far above the ring median
# Texture (smooth variation around the ring) also raises CV_excess, so
# compare samples rather than reading N_eff as an absolute grain count.
#
# Reads the raw images (never modifies them) and the outputs of main.py.
# Run from the repo root:  python -m scripts.spottiness

### Written by Mauricio Ceballos, Joester Group, Northwestern University

import numpy as np
import pandas as pd
import G2script as G2sc
from scripts.config import (DATA_DIR, PROJECT_DIR, OUTPUT_DIR, CONTROLS_FILE,
                            MASK_FILE, CTE_SETTINGS, SPOTTINESS_SETTINGS)
from scripts.analysis import RESULTS_FILE, SAMPLE_INFO_FILE

try:
    from GSASII import GSASIIimage as G2img
except ImportError:
    G2img = G2sc.G2img

SPOTTINESS_FILE = OUTPUT_DIR / "spottiness.csv"
SPOTTINESS_SUMMARY_FILE = OUTPUT_DIR / "spottiness_per_sample.csv"


def two_theta_hex(h, k, l, a, c, wavelength):
    """2theta (deg) of reflection hkl for a hexagonal cell"""
    inv_d2 = 4.0 / 3.0 * (h * h + h * k + k * k) / a**2 + l * l / c**2
    return 2 * np.degrees(np.arcsin(wavelength * np.sqrt(inv_d2) / 2))


def ring_stats(image, tth, azm, bad, tth_ring, tth_others, azm_range, settings):
    """
    Spottiness of one ring on one image.

    image: raw counts; tth, azm: 2theta and azimuth (deg) of every pixel;
    bad: True for masked pixels; tth_ring: 2theta of this reflection;
    tth_others: 2theta of other reflections (kept out of the background);
    azm_range: (start, end) of the integration range in degrees, may wrap
    past 360. Returns a dict, or None if the ring is too weak to measure.
    """
    hw = settings['ring_half_width']
    b_in, b_out = settings['bkg_inner'], settings['bkg_outer']
    start, end = azm_range
    rel = (azm - start) % 360.0                 # azimuth from the start of the range
    ok = ~bad & (rel <= end - start)

    dist = np.abs(tth - tth_ring)
    ring = ok & (dist <= hw)
    bkg = ok & (dist >= b_in) & (dist <= b_out)
    for other in tth_others:                    # keep neighbouring rings out
        bkg &= np.abs(tth - other) > b_in

    n_bins = int(np.floor((end - start) / settings['azimuth_bin']))
    edges = np.linspace(0, n_bins * settings['azimuth_bin'], n_bins + 1)
    ring_bin = np.digitize(rel[ring], edges) - 1
    bkg_bin = np.digitize(rel[bkg], edges) - 1
    ring_vals, bkg_vals = image[ring], image[bkg]

    net, var = [], []
    for b in range(n_bins):
        r = ring_vals[ring_bin == b]
        g = bkg_vals[bkg_bin == b]
        if len(r) < 5 or len(g) < 5:
            continue
        r_mean, g_mean = r.mean(), g.mean()
        net.append(r_mean - g_mean)
        # Poisson: variance of a mean of counts = mean / n
        var.append(max(r_mean, 0) / len(r) + max(g_mean, 0) / len(g))
    net, var = np.array(net), np.array(var)
    if len(net) < 10 or net.mean() <= 0:
        return None

    mean = net.mean()
    cv_obs = net.std(ddof=1) / mean
    cv_noise = np.sqrt(var.mean()) / mean
    cv_excess = np.sqrt(max(cv_obs**2 - cv_noise**2, 0.0))
    median = np.median(net)
    mad = 1.4826 * np.median(np.abs(net - median))
    spot_cut = median + settings['spot_sigma'] * max(mad, np.sqrt(var.mean()))
    return {
        'N_bins': len(net),
        'Mean_net': mean,
        'CV_obs': cv_obs,
        'CV_noise': cv_noise,
        'CV_excess': cv_excess,
        'N_eff': 1.0 / cv_excess**2 if cv_excess > 0 else np.inf,
        'Spot_fraction': float(np.mean(net > spot_cut)),
    }


def _image_arrays(image):
    """Raw counts, 2theta/azimuth maps and bad-pixel mask for a G2Image"""
    z = np.asarray(image.getImage(), dtype=float)
    controls = image.getControls()
    masks = image.getMasks()
    ta = G2img.Make2ThetaAzimuthMap(controls, (0, z.shape[0]), (0, z.shape[1]))
    position_mask = G2img.MakeMaskMap(controls, masks, (0, z.shape[0]), (0, z.shape[1]))
    if position_mask.shape != z.shape:
        raise RuntimeError(f"Mask shape {position_mask.shape} != image shape {z.shape}")
    lo, hi = masks['Thresholds'][1]
    bad = np.asarray(position_mask, dtype=bool) | ~np.isfinite(z) | (z < max(lo, 0)) | (z > hi)
    return z, ta[0], ta[1], bad, controls


def sample_spottiness(sample, scans, settings):
    """
    Spottiness of every requested scan of one sample.
    scans: rows of seq_results (T, A, C) for the scans to measure.
    """
    folder = DATA_DIR / ('AFP' if sample.startswith('AFP') else 'WT') / sample
    gpx = G2sc.G2Project(newgpx=str(PROJECT_DIR / '_spottiness.gpx'))
    geometry = None
    rows = []
    for scan in scans.itertuples():
        files = list(folder.glob(f'*_{int(scan.T)}.cbf'))
        if len(files) != 1:
            print(f"  {sample} {scan.T:g} K: expected one image, found {len(files)}")
            continue
        image = gpx.add_image(str(files[0]), fmthint='CBF')[0]
        image.loadControls(str(CONTROLS_FILE))
        image.loadMasks(str(MASK_FILE))
        controls = image.data['Image Controls']     # same tweaks as integrate_images
        controls['sag'] = 0.0
        controls['det2theta'] = 0
        if geometry is None:                        # same geometry for every scan
            z, tth, azm, bad, controls = _image_arrays(image)
            geometry = (tth, azm, bad, controls)
        else:
            z = np.asarray(image.getImage(), dtype=float)
            tth, azm, bad, controls = geometry
            bad = bad | ~np.isfinite(z) | (z < 0)
        wavelength = controls['wavelength']
        azm_range = tuple(controls['LRazimuth'])
        positions = {
            hkl: two_theta_hex(*hkl, scan.A, scan.C, wavelength)
            for hkl in settings['reflections']
        }
        for hkl, tth_ring in positions.items():
            others = [t for h, t in positions.items() if h != hkl]
            res = ring_stats(z, tth, azm, bad, tth_ring, others, azm_range, settings)
            if res is None:
                continue
            rows.append({'Sample': sample, 'T': scan.T,
                         'hkl': ''.join(str(i) for i in hkl),
                         'TwoTheta': tth_ring, **res})
    return rows


def run_spottiness(settings):
    results = pd.read_csv(RESULTS_FILE)
    info = pd.read_csv(SAMPLE_INFO_FILE)
    freeze = dict(zip(info['Sample'], info['T_freeze']))
    excluded = CTE_SETTINGS.get('excluded_scans') or {}

    rows = []
    for sample, scans in results.groupby('Sample', sort=False):
        # Same scans as the CTE fits: frozen, refined, not excluded
        scans = scans.dropna(subset=['A', 'C'])
        scans = scans[(scans['T'] < freeze.get(sample, np.inf))
                      & ~scans['T'].isin(excluded.get(sample, []))]
        print(f"Spottiness {sample}: {len(scans)} scans")
        rows += sample_spottiness(sample, scans, settings)
    per_scan = pd.DataFrame(rows)

    # One number per larva: median over scans of the mean CV_excess over
    # the measured reflections
    per_ring = per_scan.groupby(['Sample', 'T'])[['CV_excess', 'Spot_fraction']].mean()
    summary = per_ring.groupby('Sample').median().rename(columns={
        'CV_excess': 'Spottiness', 'Spot_fraction': 'Spot_fraction'}).reset_index()
    summary['N_eff'] = 1.0 / summary['Spottiness']**2
    return per_scan, summary


if __name__ == "__main__":
    per_scan, summary = run_spottiness(SPOTTINESS_SETTINGS)
    per_scan.to_csv(SPOTTINESS_FILE, index=False)
    summary.to_csv(SPOTTINESS_SUMMARY_FILE, index=False)
    pd.set_option('display.width', 200)
    print(summary.to_string(index=False))
    print(f"\nSaved {SPOTTINESS_FILE} and {SPOTTINESS_SUMMARY_FILE}")
