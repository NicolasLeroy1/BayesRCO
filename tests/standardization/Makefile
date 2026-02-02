FC = gfortran
CC = gcc
FFLAGS = -O3 -I../../new_src
CFLAGS = -O3 -I../../new_src_c -lm

FORTRAN_SRC = ../../new_src/mod_defs.f90 ../../new_src/mod_data.f90 ../../new_src/mod_io.f90 ../../new_src/mod_standardize.f90 test_standardize.f90
C_SRC = ../../new_src_c/io.c ../../new_src_c/rng.c ../../new_src_c/utils.c test_standardize.c

all: test_standardize_f test_standardize_c

test_standardize_f: $(FORTRAN_SRC)
	$(FC) $(FFLAGS) -o test_standardize_f $(FORTRAN_SRC)

test_standardize_c: $(C_SRC)
	$(CC) -O3 -I../../new_src_c -o test_standardize_c $(C_SRC) -lm

data:
	python3 generate_data.py

test: all data
	./test_standardize_f > results_f.txt
	./test_standardize_c > results_c.txt
	python3 compare_std.py results_f.txt results_c.txt

clean:
	rm -f test_standardize_f test_standardize_c results_f.txt results_c.txt *.mod test_data.*
