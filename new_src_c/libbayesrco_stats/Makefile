CC = gcc
CFLAGS = -O2 -Wall -fPIC
AR = ar
ARFLAGS = rcs

SOURCES = bayesrco_stats.c mcmc.c mcmc_additive.c mcmc_bayesCpi.c mcmc_mixture.c \
          mcmc_sampling.c mcmc_utils.c rng.c utils.c
OBJECTS = $(SOURCES:.c=.o)
TARGET = libbayesrco_stats.a

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
