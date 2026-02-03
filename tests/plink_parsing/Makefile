FC = gfortran
CC = gcc
FFLAGS = -O3 -I../../new_src
CFLAGS = -O3 -I../../new_src_c/libbayesrco_stats -I../../new_src_c/libbayesrco_io -lm

FORTRAN_SRC = ../../new_src/mod_defs.f90 ../../new_src/mod_data.f90 ../../new_src/mod_io.f90 plink_test.f90
C_SRC = plink_test.c
C_LIBS = ../../new_src_c/libbayesrco_io/libbayesrco_io.a ../../new_src_c/libbayesrco_stats/libbayesrco_stats.a

all: plink_test_f plink_test_c

plink_test_f: $(FORTRAN_SRC)
	$(FC) $(FFLAGS) -o plink_test_f $(FORTRAN_SRC)

plink_test_c: $(C_SRC)
	$(CC) $(CFLAGS) -o plink_test_c $(C_SRC) $(C_LIBS) -lm

data:
	python3 generate_data.py

test: all data
	./plink_test_f > results_f.txt
	./plink_test_c > results_c.txt
	python3 compare_plink.py results_f.txt results_c.txt

clean:
	rm -f plink_test_f plink_test_c results_f.txt results_c.txt *.mod test_data.*
