#ifndef SKY_REALISATION_H
#define SKY_REALISATION_H

#include <complex.h>

/* Default width (sigma, metres) of the Gaussian aperture footprint used to
 * turn the random uv-plane into visibilities. In wavelengths this is
 * sigma_m * frequency / c, about 1.05 at 185 MHz, which matches a primary
 * beam of FWHM ~21 degrees (MWA tile). The correlation between two
 * visibilities then falls to 50% at a separation of 2*sqrt(ln 2)*sigma,
 * about 1.75 wavelengths. */
#define SKY_SIGMA_M_DEFAULT 1.7

/* Speed of light (m/s), to turn the footprint width into wavelengths. */
#define SKY_VEL_LIGHT 299792458.0

/* Cell size (wavelengths) of the random uv-plane. */
#define SKY_CELL 0.5

/* Visibility of a simulated random sky at (u,v) in wavelengths.
 *
 * The sky is a field of independent unit-variance complex Gaussian numbers on
 * a uv grid of SKY_CELL wavelengths, Hermitian (S(-u,-v) = conj S(u,v)),
 * independent between frequency channels (chan_id), and fully determined by
 * seed. The visibility is that field seen through a Gaussian aperture
 * footprint of width sigma (wavelengths), normalised so that the visibility
 * has unit variance: <|V|^2> = 1.
 *
 * The function is pure (no state), so every observation gridded with the same
 * seed sees the same sky, and it is safe to call from several threads. */
double complex chips_sky_vis(double u, double v, long chan_id,
                             unsigned long long seed, double sigma);

/* The unit-variance complex Gaussian number of one uv cell. Exposed for tests. */
double complex chips_sky_cell(long i, long j, long chan_id,
                              unsigned long long seed);

#endif
