#!/usr/bin/env python3
"""Inspect a CMacIonize Gadget-format output: list fields, units, and
the neutral hydrogen fraction range.

Run:  python3 transient/debug_cmi_output.py <file.hdf5>
"""
import sys

import h5py
import numpy as np


def main():
    with h5py.File(sys.argv[1], 'r') as f:
        def show(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(f"{name:45s} {obj.shape} {obj.dtype}")
        f.visititems(show)
        for g in ('Units', 'Header'):
            if g in f:
                print(g, dict(f[g].attrs))
        x = f['PartType0/NeutralFractionH'][:]
        print(f"NeutralFractionH: min {x.min():.3e}, median "
              f"{np.median(x):.3e}, max {x.max():.3e}, "
              f"fraction of cells with x < 0.5: {(x < 0.5).mean():.3f}")


if __name__ == '__main__':
    main()
