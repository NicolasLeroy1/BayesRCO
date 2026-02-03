CC = gcc
CFLAGS = -O2 -Wall -fPIC -I../libbayesrco_stats
AR = ar
ARFLAGS = rcs

SOURCES = plink_io.c mcmc_core_io.c
OBJECTS = $(SOURCES:.c=.o)
TARGET = libbayesrco_io.a

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
