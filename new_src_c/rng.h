#ifndef RNG_H
#define RNG_H

#include <stdint.h>
#include <stdbool.h>

#define GFC_REAL_8_DIGITS 53

typedef struct {
    uint64_t s[4];
} prng_state;

// Core RNG functions
void manual_seed(prng_state *rs, uint64_t user_seed[4]);
void rnumber_8(double *f, uint64_t v);
uint64_t prng_next(prng_state* rs);

// Distribution functions (mirrors mod_random.f90)
double rand_uniform(prng_state *rs, double a, double b);
double rand_normal(prng_state *rs, double mean, double stdev);
double rand_exponential(prng_state *rs, double mean);
double rand_gamma(prng_state *rs, double shape, double scale);
double rand_chi_square(prng_state *rs, double dof);
double rand_scaled_inverse_chi_square(prng_state *rs, double dof, double scale);
double rand_inverse_gamma(prng_state *rs, double shape, double scale);
double rand_weibull(prng_state *rs, double shape, double scale);
double rand_cauchy(prng_state *rs, double median, double scale);
double rand_student_t(prng_state *rs, double dof);
double rand_laplace(prng_state *rs, double mean, double scale);
double rand_log_normal(prng_state *rs, double mu, double sigma);
double rand_beta(prng_state *rs, double a, double b);
void rdirichlet(prng_state *rs, int n, double *irx, double *x);
void rdirichlet2(prng_state *rs, int n, double *irx, double *x);

#endif // RNG_H
