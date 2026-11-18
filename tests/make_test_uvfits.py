"""Write a small synthetic uvfits file that gridvisdiff can read.

    python3 make_test_uvfits.py out.uvfits [n_ant] [n_time] [n_freq] [seed] [t_offset_s] [radius_m]

A random array of n_ant antennas inside a circle of radius_m (default 150 m; the
layout depends only on seed, and t_offset_s rotates it as the sky does), n_time timesteps of
2 s, n_freq channels of 80 kHz from 182 MHz, two polarisations (XX, YY),
unit-variance random visibilities, weights 1. Also writes out.uvfits.npz with
the uvw (seconds) and visibilities for checks.
"""
import sys
import numpy as np
from astropy.io import fits

out = sys.argv[1]
n_ant = int(sys.argv[2]) if len(sys.argv) > 2 else 24
n_time = int(sys.argv[3]) if len(sys.argv) > 3 else 6
n_freq = int(sys.argv[4]) if len(sys.argv) > 4 else 8
seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
t_offset = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0   # seconds, start of this observation
radius = float(sys.argv[7]) if len(sys.argv) > 7 else 150.0
rng = np.random.default_rng(seed)
C = 299792458.0
F0, DF, NPOL = 182.0e6, 80.0e3, 2
JD0 = 2460955.5

r = radius * np.sqrt(rng.random(n_ant)); th = 2 * np.pi * rng.random(n_ant)
xyz = np.stack([r * np.cos(th), r * np.sin(th), 0.3 * rng.standard_normal(n_ant)], axis=1)
a1, a2 = np.triu_indices(n_ant, 1)
nbl = len(a1)
uu, vv, ww, bl, date, vis = [], [], [], [], [], []
for t in range(n_time):
    ang = 2 * np.pi * (2.0 * t + t_offset) / 86164.0          # sky rotation in 2 s
    c, s = np.cos(ang), np.sin(ang)
    d = xyz[a2] - xyz[a1]
    u = (c * d[:, 0] - s * d[:, 1]) / C
    v = (s * d[:, 0] + c * d[:, 1]) / C
    w = d[:, 2] / C
    uu.append(u); vv.append(v); ww.append(w)
    bl.append(256.0 * (a1 + 1) + (a2 + 1))
    date.append(np.full(nbl, (0.5 + (2.0 * t + 1.0) / 86400.0)))
    vis.append((rng.standard_normal((nbl, n_freq, NPOL)) + 1j * rng.standard_normal((nbl, n_freq, NPOL))) / np.sqrt(2))
uu, vv, ww, bl, date = (np.concatenate(x).astype('f4') for x in (uu, vv, ww, bl, date))
vis = np.concatenate(vis)
data = np.zeros((len(uu), 1, 1, n_freq, NPOL, 3), 'f4')
data[:, 0, 0, :, :, 0] = vis.real; data[:, 0, 0, :, :, 1] = vis.imag; data[:, 0, 0, :, :, 2] = 1.0
g = fits.GroupData(data, parnames=['UU', 'VV', 'WW', 'BASELINE', 'DATE'], pardata=[uu, vv, ww, bl, date], bitpix=-32)
h = fits.GroupsHDU(g)
hd = h.header
hd['PZERO5'] = JD0
for i, (ct, cv, cp, cd) in enumerate([('COMPLEX', 1., 1., 1.), ('STOKES', -5., 1., -1.), ('FREQ', F0, 1., DF), ('RA', 60., 1., 1.), ('DEC', -30., 1., 1.)], start=2):
    hd[f'CTYPE{i}'] = ct; hd[f'CRVAL{i}'] = cv; hd[f'CRPIX{i}'] = cp; hd[f'CDELT{i}'] = cd
hd['OBJECT'] = 'TEST'; hd['OBSRA'] = 60.0; hd['OBSDEC'] = -30.0; hd['TELESCOP'] = 'MWA'; hd['INSTRUME'] = 'MWA'; hd['EPOCH'] = 2000.0; hd['BSCALE'] = 1.0
cols = [fits.Column('ANNAME', '8A', array=[f'T{i:03d}' for i in range(n_ant)]),
        fits.Column('STABXYZ', '3D', unit='METERS', array=xyz),
        fits.Column('NOSTA', '1J', array=np.arange(1, n_ant + 1)),
        fits.Column('MNTSTA', '1J', array=np.zeros(n_ant, 'i4')),
        fits.Column('STAXOF', '1E', unit='METERS', array=np.zeros(n_ant, 'f4')),
        fits.Column('POLTYA', '1A', array=['X'] * n_ant),
        fits.Column('POLAA', '1E', unit='DEGREES', array=np.zeros(n_ant, 'f4')),
        fits.Column('POLCALA', '1E', array=np.zeros(n_ant, 'f4')),
        fits.Column('POLTYB', '1A', array=['Y'] * n_ant),
        fits.Column('POLAB', '1E', unit='DEGREES', array=np.full(n_ant, 90, 'f4')),
        fits.Column('POLCALB', '1E', array=np.zeros(n_ant, 'f4'))]
an = fits.BinTableHDU.from_columns(cols, name='AIPS AN')
for k, v in dict(ARRAYX=-2559454.08, ARRAYY=5095372.14, ARRAYZ=-2849057.18, FREQ=F0, RDATE='2025-10-07T00:00:00', GSTIA0=0.0, DEGPDY=360.985,
                 UT1UTC=0.0, IATUTC=37.0, TIMSYS='UTC', ARRNAM='MWA', NUMORB=0, NOPCAL=0, POLARX=0.0, POLARY=0.0, FRAME='ITRF').items():
    an.header[k] = v
fits.HDUList([h, an]).writeto(out, overwrite=True)
np.savez(out + '.npz', u=uu.astype('f8'), v=vv.astype('f8'), w=ww.astype('f8'), vis=vis, n_time=n_time, nbl=nbl, f0=F0, df=DF)
print('wrote', out, 'baselines', nbl, 'times', n_time, 'chans', n_freq)
