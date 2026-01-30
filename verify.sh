#!/bin/bash
set -e

# Run Fortran version
# Run Fortran version
# Run Fortran version
echo "Running Fortran version..."
time ./build/bin/bayesRCO_fortran -bfile toy_dataset/toy_app -out toy_fortran -numit 200 -burnin 50 -thin 1 -seed 42 -catfile toy_dataset/toy_annot.txt

# Run C version
echo "Running C version..."
time ./build/bin/bayesRCO_c -bfile toy_dataset/toy_app -out toy_c -numit 200 -burnin 50 -thin 1 -seed 42 -catfile toy_dataset/toy_annot.txt

# Run R version
echo "Running R version..."
time Rscript src_R/run_bayesRCO.R -bfile toy_dataset/toy_app -out toy_r -numit 200 -seed 42 -catfile toy_dataset/toy_annot.txt

echo "All versions ran successfully!"
echo "Comparing outputs..."

# Check frequency parity (Deterministic)
if [ -f "toy_fortran.frq" ] && [ -f "toy_c.frq" ] && [ -f "toy_r.frq" ]; then
    echo "--- Frequency Comparison ---"
    diff -q toy_fortran.frq toy_c.frq && echo "Fortran vs C: Frequency files match EXACTLY." || (echo "Fortran vs C: Frequency files DIFFER:" && diff -y --suppress-common-lines toy_fortran.frq toy_c.frq | head -n 10)
    diff -q toy_fortran.frq toy_r.frq && echo "Fortran vs R: Frequency files match EXACTLY." || (echo "Fortran vs R: Frequency files DIFFER:" && diff -y --suppress-common-lines toy_fortran.frq toy_r.frq | head -n 10)
fi

# Statistical comparison of chains
if [ -f "toy_fortran.hyp" ] && [ -f "toy_c.hyp" ] && [ -f "toy_r.hyp" ]; then
    echo "--- Statistical Chain Comparison ---"
    Rscript src_R/compare_models.R toy_fortran.hyp toy_c.hyp toy_r.hyp 50
else
    echo "One or more hyperparameter chain files (.hyp) not found."
fi
