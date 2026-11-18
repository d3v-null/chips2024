#!/usr/bin/env python3
"""Reference grid for the simulated sky of `gridvisdiff -S`.

    OUTPUTDIR=/path/to/grids/ python3 sky_reference.py POL IN_EXT OUT_EXT NCHAN [SEED]

Reads   weights_POL.IN_EXT.dat   (output of prepare_diff or combine_data)
Writes  vis_tot_POL.OUT_EXT.dat, vis_diff_POL.OUT_EXT.dat and links or copies
weights_, noise_tot_ and noise_diff_ from IN_EXT to OUT_EXT.

`gridvisdiff -S` gives every visibility a simulated sky value of unit variance.
Where visibilities that see different sky modes are averaged into one cell,
part of that power is lost. This script writes what the grid would hold with
no loss at all: a value of the same variance in every cell that has data
(variance 4, because the gridder adds the two interleaved timesteps), with the
same weights. Running lssa on both and dividing,

    correction(k) = power(reference) / power(simulated sky),

gives the factor by which the power of a real sky signal is reduced by the
gridding and by the combination of observations.
"""
import os
import shutil
import sys

import numpy as np


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    pol, in_ext, out_ext, nchan = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    d = os.environ["OUTPUTDIR"]
    name = lambda kind, ext: os.path.join(d, f"{kind}_{pol}.{ext}.dat")
    if in_ext == out_ext:
        sys.exit("IN_EXT and OUT_EXT must differ")

    weights = np.memmap(name("weights", in_ext), dtype="f4", mode="r")
    if weights.size % nchan:
        sys.exit(f"weights file size {weights.size} is not a multiple of NCHAN {nchan}")
    ncell = weights.size // nchan
    rng = np.random.default_rng(seed)
    step = max(1, (1 << 24) // nchan)  # cells per chunk
    with open(name("vis_tot", out_ext), "wb") as ftot, open(name("vis_diff", out_ext), "wb") as fdiff:
        for a in range(0, ncell, step):
            b = min(ncell, a + step)
            w = np.asarray(weights[a * nchan : b * nchan])
            # complex Gaussian, <|g|^2> = 4
            g = (rng.standard_normal(w.size, dtype="f4") + 1j * rng.standard_normal(w.size, dtype="f4")) * np.float32(np.sqrt(2.0))
            g[w <= 0] = 0
            g.astype("c8").tofile(ftot)
            np.zeros(w.size, "c8").tofile(fdiff)
    for kind in ("weights", "noise_tot", "noise_diff"):
        src, dst = name(kind, in_ext), name(kind, out_ext)
        if os.path.lexists(dst):
            os.remove(dst)
        try:
            os.link(src, dst)
        except OSError:
            shutil.copyfile(src, dst)
    print(f"wrote reference grid {out_ext} ({pol}): {ncell} cells x {nchan} channels, {int((weights > 0).sum())} samples with data")


if __name__ == "__main__":
    main()
