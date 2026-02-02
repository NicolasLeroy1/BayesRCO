#include "rng.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef PI
#define PI 3.141592653589793238462
#endif

/* xoshiro256** XOR-scrambling keys for GFortran 13 */
static const uint64_t xor_keys[] = {
    0xbd0c5b6e50c2df49ULL, 0xd46061cd46e1df38ULL, 
    0xbb4f4d4ed6103544ULL, 0x114a583d0756ad39ULL
};

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

uint64_t rng_next_u64(prng_state* rs) {
    const uint64_t result = rotl(rs->s[1] * 5, 7) * 9;
    const uint64_t t = rs->s[1] << 17;

    rs->s[2] ^= rs->s[0];
    rs->s[3] ^= rs->s[1];
    rs->s[1] ^= rs->s[2];
    rs->s[0] ^= rs->s[3];

    rs->s[2] ^= t;
    rs->s[3] = rotl(rs->s[3], 45);

    return result;
}

void rng_get_double(double *f, uint64_t v) {
    uint64_t mask = ~(uint64_t)0u << (64 - GFC_REAL_8_DIGITS);
    v = v & mask;
    *f = (double)v * 5.42101086242752217e-20; 
}

void rng_seed(prng_state *rs, uint64_t user_seed[4]) {
    for (size_t i = 0; i < 4; i++)
        rs->s[i] = user_seed[i] ^ xor_keys[i];
}

// Distribution Functions

double rng_uniform(prng_state *rs, double a, double b) {
    double temp;
    rng_get_double(&temp, rng_next_u64(rs));
    return a + temp * (b - a);
}

double rng_normal(prng_state *rs, double mean, double stdev) {
    double c, r, theta, t1, t2;
    if (stdev <= 0.0) {
        return mean;
    } else {
        rng_get_double(&t1, rng_next_u64(rs));
        rng_get_double(&t2, rng_next_u64(rs));
        r = sqrt(-2.0 * log(t1));
        theta = 2.0 * PI * t2;
        c = mean + stdev * r * sin(theta);
        return c;
    }
}

double rng_exponential(prng_state *rs, double mean) {
    if (mean <= 0.0) {
        // printf("mean must be positive\n"); // Keep silent or handle error
        return 0.0;
    } else {
        double temp;
        rng_get_double(&temp, rng_next_u64(rs));
        return -mean * log(temp);
    }
}

double rng_gamma(prng_state *rs, double shape, double scale) {
    double ans, u, x, xsq, d, c, v, g, w;
    if (shape <= 0.0 || scale <= 0.0) {
        // Error handling
        return 0.0;
    }

    if (shape >= 1.0) {
        d = shape - 1.0 / 3.0;
        c = 1.0 / sqrt(9.0 * d);
        while (1) {
            x = rng_normal(rs, 0.0, 1.0);
            v = 1.0 + c * x;
            while (v <= 0.0) {
                x = rng_normal(rs, 0.0, 1.0);
                v = 1.0 + c * x;
            }
            v = v * v * v;
            rng_get_double(&u, rng_next_u64(rs));
            xsq = x * x;
            if ((u < 1.0 - .0331 * xsq * xsq) || (log(u) < 0.5 * xsq + d * (1.0 - v + log(v)))) {
                ans = scale * d * v;
                return ans;
            }
        }
    } else {
        g = rng_gamma(rs, shape + 1.0, 1.0);
        rng_get_double(&w, rng_next_u64(rs));
        ans = scale * g * pow(w, 1.0 / shape);
        return ans;
    }
}

double rng_chi_square(prng_state *rs, double dof) {
    return rng_gamma(rs, 0.5 * dof, 2.0);
}

double rng_scaled_inverse_chi_square(prng_state *rs, double dof, double scale) {
    return rng_inverse_gamma(rs, 0.5 * dof, 0.5 * dof * scale);
}

double rng_inverse_gamma(prng_state *rs, double shape, double scale) {
    return 1.0 / rng_gamma(rs, shape, 1.0 / scale);
}

double rng_weibull(prng_state *rs, double shape, double scale) {
    double temp;
    rng_get_double(&temp, rng_next_u64(rs));
    return scale * pow(-log(temp), 1.0 / shape);
}

double rng_cauchy(prng_state *rs, double median, double scale) {
    double p;
    rng_get_double(&p, rng_next_u64(rs));
    return median + scale * tan(PI * (p - 0.5));
}

double rng_student_t(prng_state *rs, double dof) {
    double y1 = rng_normal(rs, 0.0, 1.0);
    double y2 = rng_chi_square(rs, dof);
    return y1 / sqrt(y2 / dof);
}

double rng_laplace(prng_state *rs, double mean, double scale) {
    double u;
    rng_get_double(&u, rng_next_u64(rs));
    if (u < 0.5) {
        return mean + scale * log(2.0 * u);
    } else {
        return mean - scale * log(2.0 * (1.0 - u));
    }
}

double rng_log_normal(prng_state *rs, double mu, double sigma) {
    return exp(rng_normal(rs, mu, sigma));
}

double rng_beta(prng_state *rs, double a, double b) {
    double u = rng_gamma(rs, a, 1.0);
    double v = rng_gamma(rs, b, 1.0);
    return u / (u + v);
}

void rng_dirichlet(prng_state *rs, int n, double *irx, double *x) {
    double sx = 0.0;
    for (int i = 0; i < n; i++) {
        x[i] = rng_gamma(rs, irx[i], 1.0);
        sx += x[i];
    }
    for (int i = 0; i < n; i++) {
        x[i] = x[i] / sx;
    }
}

void rng_dirichlet2(prng_state *rs, int n, double *irx, double *x) {
    double sx = 0.0;
    for (int i = 0; i < n; i++) {
        if (irx[i] != 0.0) {
            x[i] = rng_gamma(rs, irx[i], 1.0);
        } else {
            x[i] = 0.0;
        }
        sx += x[i];
    }
    for (int i = 0; i < n; i++) {
        x[i] = x[i] / sx;
    }
}
