#ifndef STATS_H
#define STATS_H

#include <stddef.h>

void init_rng(int seed);
double rand_uniform(double a, double b);
double rand_normal(double mean, double stdev);
double rand_exponential(double mean);
double rand_gamma(double shape, double scale);
double rand_chi_square(double dof);
double rand_scaled_inverse_chi_square(double dof, double scale);
double rand_inverse_gamma(double shape, double scale);
double rand_beta(double a, double b);
void rdirichlet(int n, const double *alpha, double *result);

// Vector/Matrix helpers
// Vector/Matrix helpers
static inline double dot_product(const double *a, const double *b, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}
#endif // STATS_H
