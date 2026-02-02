# Makefile for BayesRCO C implementation
CC = gcc
CFLAGS = -O2 -Wall -lm

# C source files
SOURCES = main.c io.c mcmc.c mcmc_utils.c mcmc_mixture.c mcmc_additive.c mcmc_bayesCpi.c mcmc_sampling.c rng.c utils.c

# Output binary
TARGET = bayesRCO_c

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -f $(TARGET)

.PHONY: all clean
