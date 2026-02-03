FC = gfortran
CC = gcc
FFLAGS = -O3 -J.
CFLAGS = -O3 -I../../new_src_c/libbayesrco_stats -I../../new_src_c/libbayesrco_io -lm

FORTRAN_MODS = ../../new_src/mod_defs.f90 ../../new_src/mod_data.f90 ../../new_src/mod_io.f90 ../../new_src/mod_standardize.f90
C_LIBS = ../../new_src_c/libbayesrco_io/libbayesrco_io.a ../../new_src_c/libbayesrco_stats/libbayesrco_stats.a

all: test_standardize_f test_standardize_c

test_standardize_f: $(FORTRAN_MODS) test_standardize.f90
	$(FC) $(FFLAGS) -o test_standardize_f $(FORTRAN_MODS) test_standardize.f90

test_standardize_c: test_standardize.c $(C_LIBS)
	$(CC) $(CFLAGS) -o test_standardize_c test_standardize.c $(C_LIBS) -lm

$(C_LIBS):
	make -C ../../new_src_c

data:
	python3 generate_data.py

test: all data
	./test_standardize_f > results_f.txt
	./test_standardize_c > results_c.txt
	python3 compare_std.py results_f.txt results_c.txt

clean:
	rm -f test_standardize_f test_standardize_c results_f.txt results_c.txt *.mod test_data.*
