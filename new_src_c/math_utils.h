#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include "rng.h"

#define LOG_UPPER_LIMIT 700.0

/**
 * Convert log-likelihoods to normalized probabilities using log-sum-exp trick.
 * 
 * @param log_likes Input array of log-likelihoods
 * @param probs     Output array of normalized probabilities (sum to 1)
 * @param n         Number of elements
 */
void logsumexp_normalize(const double *log_likes, double *probs, int n);

/**
 * Sample an index from a discrete probability distribution.
 * 
 * @param probs Array of probabilities (should sum to 1)
 * @param n     Number of elements
 * @param rs    Random state
 * @return      Sampled index (0-based)
 */
int sample_from_probabilities(const double *probs, int n, prng_state *rs);

#endif // MATH_UTILS_H
