#include "../new_src_c/rng.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char **argv) {
    prng_state rs;
    int user_seed = 12345;
    if (argc >= 2) {
        user_seed = atoi(argv[1]);
    }

    uint32_t p[8];
    for (int i = 0; i < 8; i++) {
        p[i] = user_seed + i;
    }

    uint64_t s_packed[4];
    s_packed[0] = ((uint64_t)p[6] << 32) | (uint64_t)p[7];
    s_packed[1] = ((uint64_t)p[4] << 32) | (uint64_t)p[5];
    s_packed[2] = ((uint64_t)p[2] << 32) | (uint64_t)p[3];
    s_packed[3] = ((uint64_t)p[0] << 32) | (uint64_t)p[1];

    manual_seed(&rs, s_packed);

    printf("Uniforms:\n");
    for (int i = 0; i < 10; i++) {
        printf("%20.16f\n", rand_uniform(&rs, 0.0, 1.0));
    }

    printf("Normals:\n");
    for (int i = 0; i < 10; i++) {
        printf("%20.16f\n", rand_normal(&rs, 0.0, 1.0));
    }

    printf("Gammas (shape=2.0, scale=1.0):\n");
    for (int i = 0; i < 10; i++) {
        printf("%20.16f\n", rand_gamma(&rs, 2.0, 1.0));
    }

    return 0;
}
