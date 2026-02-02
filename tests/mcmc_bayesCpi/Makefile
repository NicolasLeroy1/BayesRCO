FC = gfortran
CC = gcc
FFLAGS = -O3 -J.
CFLAGS = -O3 -I../../new_src_c -lm

FORTRAN_MODS = ../../new_src/mod_defs.f90 ../../new_src/mod_data.f90 ../../new_src/mod_random.f90 \
               ../../new_src/mod_io.f90 ../../new_src/mod_stats.f90 ../../new_src/mod_standardize.f90 \
               ../../new_src/mod_mcmc_utils.f90 ../../new_src/mod_mcmc_bayesCpi.f90

C_SRC = ../../new_src_c/io.c ../../new_src_c/rng.c ../../new_src_c/utils.c \
        ../../new_src_c/mcmc_utils.c ../../new_src_c/mcmc_bayesCpi.c

all: test_mcmc_bayesCpi_f test_mcmc_bayesCpi_c

test_mcmc_bayesCpi_f: $(FORTRAN_MODS) test_mcmc_bayesCpi.f90
	$(FC) $(FFLAGS) -o test_mcmc_bayesCpi_f $(FORTRAN_MODS) test_mcmc_bayesCpi.f90

test_mcmc_bayesCpi_c: $(C_SRC) test_mcmc_bayesCpi.c
	$(CC) -O3 -I../../new_src_c -o test_mcmc_bayesCpi_c $(C_SRC) test_mcmc_bayesCpi.c -lm

test: all
	./test_mcmc_bayesCpi_f > results_f.txt
	./test_mcmc_bayesCpi_c > results_c.txt
	python3 compare_mcmc_bayesCpi.py results_f.txt results_c.txt

clean:
	rm -f test_mcmc_bayesCpi_f test_mcmc_bayesCpi_c results_f.txt results_c.txt *.mod
