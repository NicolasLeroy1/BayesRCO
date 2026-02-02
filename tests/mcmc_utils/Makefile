FC = gfortran
CC = gcc
FFLAGS = -O3 -J.
CFLAGS = -O3 -I../../new_src_c -lm

FORTRAN_MODS = ../../new_src/mod_defs.f90 ../../new_src/mod_data.f90 ../../new_src/mod_random.f90 \
               ../../new_src/mod_io.f90 ../../new_src/mod_stats.f90 ../../new_src/mod_standardize.f90 \
               ../../new_src/mod_mcmc_utils.f90

C_SRC = ../../new_src_c/io.c ../../new_src_c/rng.c ../../new_src_c/utils.c ../../new_src_c/mcmc_utils.c

all: test_mcmc_utils_f test_mcmc_utils_c

test_mcmc_utils_f: $(FORTRAN_MODS) test_mcmc_utils.f90
	$(FC) $(FFLAGS) -o test_mcmc_utils_f $(FORTRAN_MODS) test_mcmc_utils.f90

test_mcmc_utils_c: $(C_SRC) test_mcmc_utils.c
	$(CC) -O3 -I../../new_src_c -o test_mcmc_utils_c $(C_SRC) test_mcmc_utils.c -lm

data:
	cp ../standardization/test_data.* . 2>/dev/null || python3 ../standardization/generate_data.py

test: all data
	./test_mcmc_utils_f > results_f.txt
	./test_mcmc_utils_c > results_c.txt
	python3 compare_mcmc_utils.py results_f.txt results_c.txt

clean:
	rm -f test_mcmc_utils_f test_mcmc_utils_c results_f.txt results_c.txt *.mod test_data.*
