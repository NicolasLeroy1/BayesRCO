#!/bin/bash
set -e

# Run Fortran version
# Run Fortran version
# Run Fortran version
echo "Running Fortran version..."
time ./build/bin/bayesRCO_fortran -bfile toy_dataset/toy_app -out toy_fortran_complex -numit 50 -burnin 50 -thin 1 -seed 42 -catfile toy_dataset/complex_annot.txt

# Run C version
echo "Running C version..."
time ./build/bin/bayesRCO_c -bfile toy_dataset/toy_app -out toy_c_complex -numit 50 -burnin 50 -thin 1 -seed 42 -catfile toy_dataset/complex_annot.txt

# Run R version
echo "Running R version..."
time Rscript src_R/run_bayesRCO.R -bfile toy_dataset/toy_app -out toy_r_complex -numit 50 -seed 42 -catfile toy_dataset/complex_annot.txt

echo "All versions ran successfully!"
echo "Comparing outputs..."

# Check frequency parity (Deterministic)
# Check frequency parity (Deterministic)
if [ -f "toy_fortran_complex.frq" ] && [ -f "toy_c_complex.frq" ] && [ -f "toy_r_complex.frq" ]; then
    echo "--- Frequency Comparison ---"
    diff -q toy_fortran_complex.frq toy_c_complex.frq && echo "Fortran vs C: Frequency files match EXACTLY." || (echo "Fortran vs C: Frequency files DIFFER:" && diff -y --suppress-common-lines toy_fortran_complex.frq toy_c_complex.frq | head -n 10)
    diff -q toy_fortran_complex.frq toy_r_complex.frq && echo "Fortran vs R: Frequency files match EXACTLY." || (echo "Fortran vs R: Frequency files DIFFER:" && diff -y --suppress-common-lines toy_fortran_complex.frq toy_r_complex.frq | head -n 10)
fi

# Statistical comparison of chains
if [ -f "toy_fortran_complex.hyp" ] && [ -f "toy_c_complex.hyp" ] && [ -f "toy_r_complex.hyp" ]; then
    echo "--- Statistical Chain Comparison ---"
    Rscript src_R/compare_models.R toy_fortran_complex.hyp toy_c_complex.hyp toy_r_complex.hyp 20
else
    echo "One or more hyperparameter chain files (.hyp) not found."
fi
