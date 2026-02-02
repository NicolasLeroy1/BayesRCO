#include "math_utils.h"
#include <math.h>
#include <stdbool.h>

void logsumexp_normalize(const double *log_likes, double *probs, int n) {
    // For each element k, compute:
    // probs[k] = 1 / (1 + sum_{l != k} exp(log_likes[l] - log_likes[k]))
    // This is numerically stable due to the subtraction in the exponent.
    
    for (int k = 0; k < n; k++) {
        double skk = log_likes[k];
        double sk = 0.0;
        bool overflow = false;
        
        for (int l = 0; l < n; l++) {
            if (l == k) continue;
            double clike = log_likes[l] - skk;
            
            if (clike > LOG_UPPER_LIMIT) {
                overflow = true;
                break;
            }
            if (clike < -LOG_UPPER_LIMIT) {
                continue;
            }
            sk += exp(clike);
        }
        
        if (overflow) {
            probs[k] = 0.0;
        } else {
            probs[k] = 1.0 / (1.0 + sk);
        }
    }
}

int sample_from_probabilities(const double *probs, int n, prng_state *rs) {
    double r = rng_uniform(rs, 0.0, 1.0);
    double cumsum = 0.0;
    
    for (int k = 0; k < n; k++) {
        cumsum += probs[k];
        if (r < cumsum) {
            return k;
        }
    }
    return n - 1;  // Fallback for floating-point rounding
}
