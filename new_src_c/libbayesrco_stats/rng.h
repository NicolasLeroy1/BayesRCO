#ifndef RNG_H
#define RNG_H

#include <stdint.h>
#include <stdbool.h>

#define GFC_REAL_8_DIGITS 53

typedef struct {
    uint64_t s[4];
} prng_state;

// Core RNG functions
void rng_seed(prng_state *rs, uint64_t user_seed[4]);
void rng_get_double(double *f, uint64_t v);
uint64_t rng_next_u64(prng_state* rs);

// Distribution functions (mirrors mod_random.f90)
double rng_uniform(prng_state *rs, double a, double b);
double rng_normal(prng_state *rs, double mean, double stdev);
double rng_exponential(prng_state *rs, double mean);
double rng_gamma(prng_state *rs, double shape, double scale);
double rng_chi_square(prng_state *rs, double dof);
double rng_scaled_inverse_chi_square(prng_state *rs, double dof, double scale);
double rng_inverse_gamma(prng_state *rs, double shape, double scale);
double rng_weibull(prng_state *rs, double shape, double scale);
double rng_cauchy(prng_state *rs, double median, double scale);
double rng_student_t(prng_state *rs, double dof);
double rng_laplace(prng_state *rs, double mean, double scale);
double rng_log_normal(prng_state *rs, double mu, double sigma);
double rng_beta(prng_state *rs, double a, double b);
void rng_dirichlet(prng_state *rs, int n, double *irx, double *x);
void rng_dirichlet2(prng_state *rs, int n, double *irx, double *x);

/* Backward-compatible aliases for legacy code */
#define rand_uniform rng_uniform
#define rand_normal rng_normal
#define rand_gamma rng_gamma
#define rand_chi_square rng_chi_square
#define rand_scaled_inverse_chi_square rng_scaled_inverse_chi_square
#define rdirichlet rng_dirichlet
#define rdirichlet2 rng_dirichlet2

#endif // RNG_H
