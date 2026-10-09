"""Download and clean public ZTF DR light curves for the transient quasar.

Queries the IRSA ZTF light-curve API (no account needed) at the object
position, cleans each band, converts AB magnitudes to flux in mJy, and
saves per-band npz files under transient/data/.

Run: python3 transient/fetch_ztf_lc.py
"""

import io
import os

import numpy as np

RA = 25.1172
DEC = 13.3891

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'data')

API_URL = 'https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves'
RADIUS_DEG = 0.000833  # 3 arcsec
BANDS = {'zg': 'g', 'zr': 'r', 'zi': 'i'}
MAGERR_MAX = 0.2


def fetch_raw_table():
    import pandas as pd
    import requests

    params = {
        'POS': f'CIRCLE {RA} {DEC} {RADIUS_DEG}',
        'BANDNAME': 'g,r,i',
        'FORMAT': 'csv',
    }
    print(f'Querying IRSA ZTF light-curve API at RA={RA}, Dec={DEC} ...')
    resp = requests.get(API_URL, params=params, timeout=300)
    resp.raise_for_status()
    df = pd.read_csv(io.StringIO(resp.text))
    os.makedirs(OUTDIR, exist_ok=True)
    raw_path = os.path.join(OUTDIR, 'ztf_irsa_raw.csv')
    df.to_csv(raw_path, index=False)
    print(f'saved {raw_path} ({len(df)} rows, columns: {list(df.columns)})')
    return df


def clean_and_save(df):
    for fcode, band in BANDS.items():
        sub = df[df['filtercode'] == fcode]
        n_all = len(sub)
        sub = sub[(sub['catflags'] == 0)
                  & np.isfinite(sub['mag']) & np.isfinite(sub['magerr'])
                  & (sub['magerr'] > 0) & (sub['magerr'] < MAGERR_MAX)]
        if len(sub) == 0:
            print(f'{band}: no points survive cleaning ({n_all} raw), skipped')
            continue
        nfields = sub['oid'].nunique() if 'oid' in sub.columns else 1
        order = np.argsort(sub['mjd'].values)
        mjd = sub['mjd'].values[order]
        mag = sub['mag'].values[order]
        magerr = sub['magerr'].values[order]
        flux_mjy = 3.631e6 * 10.0 ** (-0.4 * mag)
        fluxerr_mjy = flux_mjy * 0.4 * np.log(10.0) * magerr
        path = os.path.join(OUTDIR, f'lc_ztf_{band}.npz')
        np.savez(path, mjd=mjd, mag=mag, magerr=magerr,
                 flux_mjy=flux_mjy, fluxerr_mjy=fluxerr_mjy, band=band)
        print(f'saved {path} ({len(mjd)} of {n_all} points, {nfields} object IDs, '
              f'MJD {mjd.min():.1f}-{mjd.max():.1f}, '
              f'<mag> = {mag.mean():.2f})')


if __name__ == '__main__':
    clean_and_save(fetch_raw_table())
