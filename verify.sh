#!/bin/bash
set -e

# Create temporary directory
TMP_DIR="temp_verify"
mkdir -p $TMP_DIR
echo "Using temporary directory: $TMP_DIR"

# Cleanup function
cleanup() {
    echo "Cleaning up temporary directory..."
    rm -rf $TMP_DIR
}
trap cleanup EXIT

# Run Fortran version
echo "Running Fortran version..."
time ./build/bin/bayesRCO_fortran -bfile toy_dataset/toy_app -out $TMP_DIR/toy_fortran -numit 200 -burnin 0 -thin 1 -seed 42 -catfile toy_dataset/toy_annot.txt

# Run Old Fortran version
echo "Running Old Fortran version..."
time ./build/bin/bayesRCO_old_fortran -bfile toy_dataset/toy_app -out $TMP_DIR/toy_old_fortran -numit 200 -burnin 0 -thin 1 -seed 42 -catfile toy_dataset/toy_annot.txt

# Run C version
echo "Running C version..."
time ./build/bin/bayesRCO_c -bfile toy_dataset/toy_app -out $TMP_DIR/toy_c -numit 200 -burnin 0 -thin 1 -seed 42 -catfile toy_dataset/toy_annot.txt

# Run R version
echo "Running R version..."
time Rscript src_R/run_bayesRCO.R -bfile toy_dataset/toy_app -out $TMP_DIR/toy_r -numit 200 -burnin 0 -seed 42 -catfile toy_dataset/toy_annot.txt

# Run Deduplicated Old Fortran version
echo "Running Deduplicated Old Fortran version..."
time ./build/bin/bayesRCO_deduplicated -bfile toy_dataset/toy_app -out $TMP_DIR/toy_dedup -numit 200 -burnin 0 -thin 1 -seed 42 -catfile toy_dataset/toy_annot.txt

# Run R Mimic version
echo "Running R Mimic version..."
time Rscript src_R_mimic/run_bayesR_mimic.R -bfile toy_dataset/toy_app -out $TMP_DIR/toy_r_mimic -numit 200 -burnin 0 -seed 42 -catfile toy_dataset/toy_annot.txt

echo "All versions ran successfully!"
echo "Comparing outputs..."

# Check frequency parity (Deterministic)
if [ -f "$TMP_DIR/toy_fortran.frq" ] && [ -f "$TMP_DIR/toy_c.frq" ] && [ -f "$TMP_DIR/toy_r.frq" ] && [ -f "$TMP_DIR/toy_old_fortran.frq" ]; then
    echo "--- Frequency Comparison ---"
    diff -q $TMP_DIR/toy_fortran.frq $TMP_DIR/toy_c.frq && echo "Fortran vs C: Frequency files match EXACTLY." || (echo "Fortran vs C: Frequency files DIFFER:" && diff -y --suppress-common-lines $TMP_DIR/toy_fortran.frq $TMP_DIR/toy_c.frq | head -n 10)
    diff -q $TMP_DIR/toy_fortran.frq $TMP_DIR/toy_r.frq && echo "Fortran vs R: Frequency files match EXACTLY." || (echo "Fortran vs R: Frequency files DIFFER:" && diff -y --suppress-common-lines $TMP_DIR/toy_fortran.frq $TMP_DIR/toy_r.frq | head -n 10)
    diff -q $TMP_DIR/toy_fortran.frq $TMP_DIR/toy_old_fortran.frq && echo "Fortran vs Old Fortran: Frequency files match EXACTLY." || (echo "Fortran vs Old Fortran: Frequency files DIFFER:" && diff -y --suppress-common-lines $TMP_DIR/toy_fortran.frq $TMP_DIR/toy_old_fortran.frq | head -n 10)
fi

# Statistical comparison of chains
if [ -f "$TMP_DIR/toy_fortran.hyp" ] && [ -f "$TMP_DIR/toy_c.hyp" ] && [ -f "$TMP_DIR/toy_r.hyp" ] && [ -f "$TMP_DIR/toy_old_fortran.hyp" ] && [ -f "$TMP_DIR/toy_r_mimic.hyp" ]; then
    echo "--- Statistical Chain Comparison ---"
    # Arguments: <C_hyp> <Fortran_New_hyp> <Fortran_Old_hyp> <Fortran_Dedup_hyp> <R_hyp> <R_Mimic_hyp> <Burnin> <OutputPDF>
    Rscript src_R/compare_models.R $TMP_DIR/toy_fortran.hyp $TMP_DIR/toy_c.hyp $TMP_DIR/toy_r.hyp $TMP_DIR/toy_old_fortran.hyp $TMP_DIR/toy_dedup.hyp $TMP_DIR/toy_r_mimic.hyp 50 similarity_visualization.pdf

else
    echo "One or more hyperparameter chain files (.hyp) not found."
fi

