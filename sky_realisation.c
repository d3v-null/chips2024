/* A reproducible random sky, sampled at arbitrary (u,v).
 *
 * Used by gridvisdiff -S to grid a simulated sky signal through the real
 * sampling, weights and flags of an observation, so that the fraction of sky
 * power that survives gridding and the combination of observations can be
 * measured instead of assumed. See README.md, "Measuring the signal
 * normalisation".
 */
#include <math.h>
#include <stdint.h>
#include "sky_realisation.h"

/* splitmix64 finaliser: a well-mixed 64-bit hash */
static inline uint64_t mix64(uint64_t z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/* uniform in (0,1], 53 bits */
static inline double to_unit(uint64_t x) {
    return ((double)(x >> 11) + 1.0) * (1.0 / 9007199254740992.0);
}

double complex chips_sky_cell(long i, long j, long chan_id,
                              unsigned long long seed) {
    int conjugate = 0;
    uint64_t h, a, b;
    double r, phi;

    /* Hermitian: the cell at (-i,-j) is the conjugate of the cell at (i,j).
     * Canonical half: i > 0, or i == 0 and j >= 0. */
    if ((i < 0) || ((i == 0) && (j < 0))) {
        i = -i;
        j = -j;
        conjugate = 1;
    }
    h = mix64((uint64_t)seed);
    h = mix64(h ^ (uint64_t)(int64_t)chan_id);
    h = mix64(h ^ (uint64_t)(int64_t)i);
    h = mix64(h ^ (uint64_t)(int64_t)j);
    a = mix64(h);
    b = mix64(a);

    /* Box-Muller: complex Gaussian with <|g|^2> = 1 */
    r = sqrt(-log(to_unit(a)));
    phi = 2.0 * M_PI * to_unit(b);
    if ((i == 0) && (j == 0)) {
        /* the origin is its own conjugate: real, unit variance */
        return r * M_SQRT2 * cos(phi);
    }
    return conjugate ? r * cos(phi) - I * r * sin(phi)
                     : r * cos(phi) + I * r * sin(phi);
}

#define SKY_MAX_HALF 64

double complex chips_sky_vis(double u, double v, long chan_id,
                             unsigned long long seed, double sigma) {
    long i0 = lround(u / SKY_CELL), j0 = lround(v / SKY_CELL);
    int half = (int)ceil(3.0 * sigma / SKY_CELL), a, b;
    double wu[2 * SKY_MAX_HALF + 1], wv[2 * SKY_MAX_HALF + 1];
    double su = 0.0, sv = 0.0, d;
    double complex sum = 0.0;

    if (half < 1) half = 1;
    if (half > SKY_MAX_HALF) half = SKY_MAX_HALF;

    /* separable Gaussian footprint centred on (u,v) */
    for (a = -half; a <= half; a++) {
        d = ((double)(i0 + a)) * SKY_CELL - u;
        wu[a + half] = exp(-0.5 * d * d / (sigma * sigma));
        su += wu[a + half] * wu[a + half];
        d = ((double)(j0 + a)) * SKY_CELL - v;
        wv[a + half] = exp(-0.5 * d * d / (sigma * sigma));
        sv += wv[a + half] * wv[a + half];
    }
    for (a = -half; a <= half; a++) {
        double complex row = 0.0;
        for (b = -half; b <= half; b++) {
            row += wv[b + half] * chips_sky_cell(i0 + a, j0 + b, chan_id, seed);
        }
        sum += wu[a + half] * row;
    }
    /* unit variance: divide by the root of the sum of squared weights */
    return sum / sqrt(su * sv);
}
