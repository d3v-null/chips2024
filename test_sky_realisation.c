/* Checks of the simulated sky used by gridvisdiff -S.
 *
 *   make test_sky && ./test_sky_realisation
 *
 * Exit status 0 if every check passes. With an argument, also prints the
 * measured correlation against separation as CSV on stdout.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "sky_realisation.h"

static int failures = 0;

static void check(const char *what, double got, double want, double tol) {
    int ok = fabs(got - want) <= tol;
    fprintf(stderr, "%s %-58s got %9.5f  want %9.5f  (tol %.3f)\n",
            ok ? "ok  " : "FAIL", what, got, want, tol);
    if (!ok) failures++;
}

int main(int argc, char *argv[]) {
    const unsigned long long seed = 20251007ULL;
    const double sigma = 1.05; /* wavelengths */
    const int n = 40000;
    double seps[] = {0.0, 0.25, 0.5, 1.0, 1.75, 2.5, 3.5, 5.0, 8.0};
    int nsep = sizeof(seps) / sizeof(seps[0]), s, k;
    double var = 0.0, re2 = 0.0, im2 = 0.0, reim = 0.0;
    double complex pseudo = 0.0;

    /* points spread over 10..300 wavelengths, any direction, 64 channels */
    srand48(1);
    if (argc > 1) printf("separation_wavelengths,correlation_measured,correlation_expected\n");
    for (s = 0; s < nsep; s++) {
        double complex c12 = 0.0;
        double p1 = 0.0, p2 = 0.0, d = seps[s];
        for (k = 0; k < n; k++) {
            double r = 10.0 + 290.0 * drand48(), th = 2.0 * M_PI * drand48(), ph = 2.0 * M_PI * drand48();
            double u = r * cos(th), v = r * sin(th);
            long ch = 2000 + (long)(64.0 * drand48());
            double complex a = chips_sky_vis(u, v, ch, seed, sigma);
            double complex b = chips_sky_vis(u + d * cos(ph), v + d * sin(ph), ch, seed, sigma);
            c12 += a * conj(b);
            p1 += creal(a * conj(a));
            p2 += creal(b * conj(b));
            if (s == 0) {
                var += creal(a * conj(a));
                re2 += creal(a) * creal(a);
                im2 += cimag(a) * cimag(a);
                reim += creal(a) * cimag(a);
                pseudo += a * a;
            }
        }
        {
            char label[128];
            double got = creal(c12) / sqrt(p1 * p2), want = exp(-d * d / (4.0 * sigma * sigma));
            snprintf(label, sizeof(label), "correlation at separation %.2f wavelengths", d);
            check(label, got, want, 0.015);
            if (argc > 1) printf("%.2f,%.5f,%.5f\n", d, got, want);
        }
    }
    check("variance <|V|^2>", var / n, 1.0, 0.01);
    check("variance of the real part", re2 / n, 0.5, 0.01);
    check("variance of the imaginary part", im2 / n, 0.5, 0.01);
    check("real-imaginary covariance", reim / n, 0.0, 0.01);
    check("pseudo-variance |<V V>|", cabs(pseudo) / n, 0.0, 0.01);

    /* Hermitian, deterministic, seed and channel dependent */
    {
        double herm = 0.0, same = 0.0;
        double complex cc = 0.0, cs = 0.0;
        double pc1 = 0.0, pc2 = 0.0, ps2 = 0.0;
        for (k = 0; k < 20000; k++) {
            double u = 300.0 * (drand48() - 0.5), v = 300.0 * (drand48() - 0.5);
            double complex a = chips_sky_vis(u, v, 2100, seed, sigma);
            double complex b = chips_sky_vis(-u, -v, 2100, seed, sigma);
            double complex c = chips_sky_vis(u, v, 2101, seed, sigma);
            double complex e = chips_sky_vis(u, v, 2100, seed + 1, sigma);
            herm += cabs(a - conj(b));
            same += cabs(a - chips_sky_vis(u, v, 2100, seed, sigma));
            cc += a * conj(c);
            cs += a * conj(e);
            pc1 += creal(a * conj(a));
            pc2 += creal(c * conj(c));
            ps2 += creal(e * conj(e));
        }
        check("Hermitian: mean |V(u,v) - conj V(-u,-v)|", herm / 20000, 0.0, 1e-12);
        check("deterministic: mean |V - V again|", same / 20000, 0.0, 0.0);
        check("correlation between neighbouring channels", cabs(cc) / sqrt(pc1 * pc2), 0.0, 0.03);
        check("correlation between two seeds", cabs(cs) / sqrt(pc1 * ps2), 0.0, 0.03);
    }
    /* same answer whichever cell the rounding picks: continuity across a cell edge */
    {
        double jump = 0.0;
        for (k = 0; k < 20000; k++) {
            double u = SKY_CELL * (floor(400.0 * drand48()) - 200.0 + 0.5), v = 300.0 * (drand48() - 0.5);
            jump += cabs(chips_sky_vis(u - 1e-9, v, 2100, seed, sigma) - chips_sky_vis(u + 1e-9, v, 2100, seed, sigma));
        }
        check("continuity across a cell edge: mean jump", jump / 20000, 0.0, 0.02);
    }
    fprintf(stderr, failures ? "%d check(s) FAILED\n" : "all checks passed\n", failures);
    return failures ? 1 : 0;
}
