"""Check the grid written by gridvisdiff -S against the simulated sky evaluated directly.

    cc -O2 -shared -fPIC -o libsky.so ../sky_realisation.c -lm
    python3 make_test_uvfits.py t.uvfits
    OUTPUTDIR=$PWD/out/ ../gridvisdiff t.uvfits 1000000001 sky 1 -f 1 -u 150 -S 7
    python3 check_sky_grid.py out sky t.uvfits 7

For baselines with no neighbour within 12 wavelengths, the cell under the
baseline must hold the simulated sky at that (u,v), including baselines that
the gridder flips to the conjugate position (w < 0).
"""
import sys, ctypes, numpy as np
S = sys.argv[1]; ext = sys.argv[2]; uvf = sys.argv[3]; seed = int(sys.argv[4]); umax = 150.
lib = ctypes.CDLL('./libsky.so')
class C(ctypes.Structure): _fields_ = [('re', ctypes.c_double), ('im', ctypes.c_double)]
lib.chips_sky_vis.restype = C
lib.chips_sky_vis.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_long, ctypes.c_ulonglong, ctypes.c_double]
def sky(u, v, ch, sig):
    r = lib.chips_sky_vis(u, v, ch, seed, sig); return complex(r.re, r.im)
d = np.load(uvf + '.npz'); nbl = int(d['nbl']); nt = int(d['n_time']); f0 = float(d['f0']); df = float(d['df'])
us = 2 * int(umax / 0.5); nf = d['vis'].shape[1]
for pol in ('xx', 'yy'):
    bv = np.fromfile(f'{S}/bv_{pol}.{ext}.dat', 'c16').reshape(us * us, nf)
    w = np.fromfile(f'{S}/weightsc_{pol}.{ext}.dat', 'c16').reshape(us * us, nf).real
    G = np.where(w > 0, bv / (2 * np.where(w > 0, w, 1)), 0)
    print(pol, 'cells with data', int((w[:, 0] > 0).sum()), 'mean |G|^2 over cells with data', float((abs(G[w > 0]) ** 2).mean()))
    # isolated baselines: nearest other baseline (either sign) further than 12 wavelengths at the top channel
    u0 = d['u'][:nbl]; v0 = d['v'][:nbl]; w0 = d['w'][:nbl]
    errs = []; nflip = 0; used = 0
    for ch in (0, nf - 1):
        f = f0 + ch * df; sig = 1.7 * f / 299792458.0; chan = int(round(f / df))
        U = u0 * f; V = v0 * f
        P = np.concatenate([np.stack([U, V], 1), -np.stack([U, V], 1)])
        for b in range(nbl):
            dist = np.hypot(P[:, 0] - U[b], P[:, 1] - V[b]); dist[b] = 1e9
            if dist.min() < 12 or np.hypot(U[b], V[b]) > 0.9 * umax: continue
            flip = w0[b] < 0
            ug, vg = (-U[b], -V[b]) if flip else (U[b], V[b])
            ii = int(round(ug / 0.5)) + us // 2 - 1; jj = int(round(vg / 0.5)) + us // 2 - 1
            g = G[ii * us + jj, ch]
            # mean over the timesteps of the simulated sky at the gridded position
            exp = np.mean([sky((-1 if flip else 1) * d['u'][t * nbl + b] * f, (-1 if flip else 1) * d['v'][t * nbl + b] * f, chan, sig) for t in range(nt)])
            errs.append(abs(g - exp)); nflip += flip; used += 1
    errs = np.array(errs)
    print(pol, 'isolated baselines checked', used, 'of which flipped (w<0)', int(nflip), 'max |G - expected|', errs.max(), 'median', np.median(errs))
