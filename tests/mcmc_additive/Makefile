FC = gfortran
CC = gcc
FFLAGS = -O3 -J.
CFLAGS = -O3 -I../../new_src_c/libbayesrco_stats -I../../new_src_c/libbayesrco_io -lm

FORTRAN_MODS = ../../new_src/mod_defs.f90 ../../new_src/mod_data.f90 ../../new_src/mod_random.f90 \
               ../../new_src/mod_io.f90 ../../new_src/mod_stats.f90 ../../new_src/mod_standardize.f90 \
               ../../new_src/mod_mcmc_utils.f90 ../../new_src/mod_mcmc_additive.f90

C_LIBS = ../../new_src_c/libbayesrco_io/libbayesrco_io.a ../../new_src_c/libbayesrco_stats/libbayesrco_stats.a

all: test_mcmc_additive_f test_mcmc_additive_c

test_mcmc_additive_f: $(FORTRAN_MODS) test_mcmc_additive.f90
	$(FC) $(FFLAGS) -o test_mcmc_additive_f $(FORTRAN_MODS) test_mcmc_additive.f90

test_mcmc_additive_c: test_mcmc_additive.c $(C_LIBS)
	$(CC) $(CFLAGS) -o test_mcmc_additive_c test_mcmc_additive.c $(C_LIBS) -lm

$(C_LIBS):
	make -C ../../new_src_c

data:
	cp ../standardization/test_data.* . 2>/dev/null || python3 ../standardization/generate_data.py

test: all data
	./test_mcmc_additive_f > results_f.txt
	./test_mcmc_additive_c > results_c.txt
	python3 compare_mcmc_additive.py results_f.txt results_c.txt

clean:
	rm -f test_mcmc_additive_f test_mcmc_additive_c results_f.txt results_c.txt *.mod test_data.*
