#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "stats.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void init_rng(int seed) {
    if (seed == 0) {
        srand(time(NULL));
    } else {
        srand(seed);
    }
}

static double rng_double() {
    return (double)rand() / ((double)RAND_MAX + 1.0);
}

double rand_uniform(double a, double b) {
    return a + rng_double() * (b - a);
}

double rand_normal(double mean, double stdev) {
    if (stdev <= 0.0) return mean;
    
    static double z1;
    static int generate = 0;
    generate = !generate;

    if (!generate) return mean + stdev * z1;

    double u1, u2;
    do {
        u1 = rng_double();
        u2 = rng_double();
    } while (u1 <= 1e-7);

    double r = sqrt(-2.0 * log(u1));
    double theta = 2.0 * M_PI * u2;
    z1 = r * sin(theta);
    return mean + stdev * r * cos(theta);
}

double rand_exponential(double mean) {
    double u;
    do { u = rng_double(); } while (u <= 1e-7);
    return -mean * log(u);
}

double rand_gamma(double shape, double scale) {
    if (shape < 1.0) {
        return rand_gamma(shape + 1.0, scale) * pow(rng_double(), 1.0 / shape);
    }

    double d = shape - 1.0 / 3.0;
    double c = 1.0 / sqrt(9.0 * d);
    double x, v, u;

    while (1) {
        do {
            x = rand_normal(0.0, 1.0);
            v = 1.0 + c * x;
        } while (v <= 0.0);

        v = v * v * v;
        u = rng_double();

        if (u < 1.0 - 0.0331 * pow(x, 4)) return scale * d * v;
        if (log(u) < 0.5 * x * x + d * (1.0 - v + log(v))) return scale * d * v;
    }
}

double rand_chi_square(double dof) {
    return rand_gamma(0.5 * dof, 2.0);
}

double rand_inverse_gamma(double shape, double scale) {
    return 1.0 / rand_gamma(shape, 1.0 / scale);
}

double rand_scaled_inverse_chi_square(double dof, double scale) {
    return rand_inverse_gamma(0.5 * dof, 0.5 * dof * scale);
}

double rand_beta(double a, double b) {
    double u = rand_gamma(a, 1.0);
    double v = rand_gamma(b, 1.0);
    return u / (u + v);
}

void rdirichlet(int n, const double *alpha, double *result) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        result[i] = rand_gamma(alpha[i], 1.0);
        sum += result[i];
    }
    for (int i = 0; i < n; i++) {
        result[i] /= sum;
    }
}


