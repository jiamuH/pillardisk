#!/usr/bin/env python3
"""
fetch_spectra.py - download the public spectra of the transient quasar
(Liu et al. "Main Target": R.A. = 25.1172, Dec. = +13.3891, z = 0.494).

Epochs:
  - SDSS Legacy   2000-12-04 (MJD ~ 51876)   via astroquery.sdss
  - SDSS/eBOSS    2019-02-09 (MJD ~ 58523)   via astroquery.sdss
  - DESI DR1      2021-10-31 / 11-02         via SPARCL (sparclclient)

Each epoch is saved as transient/data/epoch_<name>.npz with fields:
  wave_obs [Angstrom, observed frame], flam [1e-17 erg/s/cm2/A, observed
  frame], ivar [inverse variance of flam], mjd, label.

Run:  python3 transient/fetch_spectra.py
"""

import os

import numpy as np

RA = 25.1172
DEC = 13.3891
REDSHIFT = 0.494
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


def save_epoch(name, wave, flam, ivar, mjd, label):
    os.makedirs(OUTDIR, exist_ok=True)
    good = np.isfinite(wave) & np.isfinite(flam) & (wave > 0)
    path = os.path.join(OUTDIR, f'epoch_{name}.npz')
    np.savez(path, wave_obs=wave[good], flam=flam[good],
             ivar=np.where(np.isfinite(ivar[good]), ivar[good], 0.0),
             mjd=mjd, label=label)
    print(f"  saved {path}  ({good.sum()} pixels, "
          f"{wave[good].min():.0f}-{wave[good].max():.0f} A, MJD {mjd})")


def fetch_sdss():
    """Both SDSS epochs (Legacy 2000 + eBOSS 2019) by cone search."""
    from astroquery.sdss import SDSS
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    co = SkyCoord(RA, DEC, unit='deg')
    xid = None
    for dr in (17, 16):
        try:
            xid = SDSS.query_region(co, radius=5 * u.arcsec, spectro=True,
                                    data_release=dr)
        except Exception as exc:
            print(f"  SDSS query (DR{dr}) failed: {exc}")
            continue
        if xid is not None and len(xid) > 0:
            print(f"  DR{dr}: found {len(xid)} spectroscopic match(es)")
            break
    if xid is None or len(xid) == 0:
        print("  no SDSS spectra found; skipping SDSS epochs")
        return

    seen = set()
    for row in xid:
        key = (int(row['plate']), int(row['mjd']), int(row['fiberID']))
        if key in seen:
            continue
        seen.add(key)
        plate, mjd, fiber = key
        try:
            hdus = SDSS.get_spectra(plate=plate, mjd=mjd, fiberID=fiber,
                                    data_release=17)[0]
        except Exception as exc:
            print(f"  download failed for plate {plate} mjd {mjd} "
                  f"fiber {fiber}: {exc}")
            continue
        coadd = hdus[1].data
        wave = 10.0 ** coadd['loglam']
        name = 'sdss2000' if mjd < 55000 else 'eboss2019'
        label = 'SDSS-2000' if mjd < 55000 else 'SDSS-2019'
        save_epoch(name, wave, coadd['flux'], coadd['ivar'], mjd, label)


def fetch_desi():
    """DESI DR1 spectrum via SPARCL. Uses the coadded (brz-combined) flux."""
    try:
        from sparcl.client import SparclClient
    except ImportError:
        print("  sparclclient not installed; skipping DESI epoch "
              "(pip install sparclclient)")
        return
    try:
        client = SparclClient()
        found = client.find(
            outfields=['sparcl_id', 'ra', 'dec', 'redshift', 'data_release',
                       'spectype'],
            constraints={'data_release': ['DESI-DR1'],
                         'ra': [RA - 0.002, RA + 0.002],
                         'dec': [DEC - 0.002, DEC + 0.002]})
    except Exception as exc:
        print(f"  SPARCL query failed: {exc}; skipping DESI epoch")
        return
    if len(found.records) == 0:
        print("  no DESI DR1 record found; skipping DESI epoch")
        return
    ids = [rec['sparcl_id'] for rec in found.records]
    print(f"  DESI DR1: found {len(ids)} record(s)")
    try:
        res = client.retrieve(uuid_list=ids,
                              include=['sparcl_id', 'wavelength', 'flux',
                                       'ivar', 'redshift'])
    except Exception as exc:
        print(f"  SPARCL retrieve failed: {exc}; skipping DESI epoch")
        return
    for k, rec in enumerate(res.records):
        name = 'desi2021' if k == 0 else f'desi2021_{k}'
        save_epoch(name, np.asarray(rec.wavelength), np.asarray(rec.flux),
                   np.asarray(rec.ivar), 59518, 'DESI-2021')


if __name__ == '__main__':
    print("Fetching SDSS epochs (Legacy 2000 + eBOSS 2019)...")
    fetch_sdss()
    print("Fetching DESI DR1 epoch...")
    fetch_desi()
    print("Done. The figure driver uses whatever epochs are present in "
          "transient/data/.")
