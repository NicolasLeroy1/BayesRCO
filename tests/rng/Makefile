FC = gfortran
CC = gcc
FFLAGS = -O3 -I../../new_src
CFLAGS = -O3 -I../../new_src_c/libbayesrco_stats -lm

FORTRAN_SRC = ../../new_src/mod_defs.f90 ../../new_src/mod_random.f90 rng_test.f90
C_SRC = rng_test.c
C_LIBS = ../../new_src_c/libbayesrco_stats/libbayesrco_stats.a

all: rng_test_f rng_test_c

rng_test_f: $(FORTRAN_SRC)
	$(FC) $(FFLAGS) -o rng_test_f $(FORTRAN_SRC)

rng_test_c: $(C_SRC)
	$(CC) $(CFLAGS) -o rng_test_c $(C_SRC) $(C_LIBS) -lm

test: all
	./rng_test_f 12345 > results_f.txt
	./rng_test_c 12345 > results_c.txt
	python3 compare_rng.py results_f.txt results_c.txt

clean:
	rm -f rng_test_f rng_test_c results_f.txt results_c.txt *.mod
