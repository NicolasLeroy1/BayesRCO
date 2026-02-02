#!/bin/bash

# Configuration
BIN_DIR="bin"
TOY_DIR="toy_dataset"
NUMIT=300
THIN=5
BURNIN=0
SEED=10
NCAT=4

# Binaries
BIN_REF="$BIN_DIR/bayesRCO_ref"
BIN_NEW="$BIN_DIR/bayesRCO"
BIN_C="$BIN_DIR/bayesRCO_c"

# Output directories (organized under outputs/)
OUT_REF="outputs/ref"
OUT_NEW="outputs/new"
OUT_C="outputs/c"

# Root directory for absolute paths
ROOT_DIR=$(pwd)

# Functions
run_model() {
    local bin=$1
    local out_dir=$2
    local label=$3
    local prefix=$4
    
    echo "--- Running $label ---"
    mkdir -p "$out_dir"
    cp $TOY_DIR/toy_annot.txt $TOY_DIR/toy_app.* $TOY_DIR/toy_test.* "$out_dir/"
    
    (
        cd "$out_dir" || exit 1
        # Run training using absolute path for binary
        "$ROOT_DIR/$bin" -bfile toy_app -out "${prefix}_train" -seed $SEED -ncat $NCAT -catfile toy_annot.txt -burnin $BURNIN -numit $NUMIT -thin $THIN
    )
}

compare_exact() {
    local dir1=$1
    local dir2=$2
    local prefix1=$3
    local prefix2=$4
    local label=$5
    local diff_found=0

    echo "--- Comparing Exact Parity: $label ($dir1 vs $dir2) ---"
    for ext in model param frq hyp; do 
        file1="$dir1/${prefix1}_train.$ext"
        file2="$dir2/${prefix2}_train.$ext"

        if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
            echo "Error: Missing files for .$ext comparison"
            diff_found=1
            continue
        fi

        if diff -q "$file1" "$file2" > /dev/null; then
            echo "[PASS] .$ext files are identical"
        else
            echo "[FAIL] .$ext files differ"
            diff "$file1" "$file2" | head -n 5
            diff_found=1
        fi
    done
    return $diff_found
}

# Main Execution
echo "Starting BayesRCO Verification System (REF vs NEW vs C)"

# Build using centralized Makefile
echo "Building versions..."
make clean
make all || { echo "Build failed"; exit 1; }

# Run versions
run_model "$BIN_REF" "$OUT_REF" "Reference (src/)" "ref"
run_model "$BIN_NEW" "$OUT_NEW" "New Modular (new_src/)" "new"
run_model "$BIN_C" "$OUT_C" "C Port (new_src_c/)" "c"

# Compare
EXIT_CODE=0

# Exact comparison
echo ""
compare_exact "$OUT_REF" "$OUT_NEW" "ref" "new" "REF vs NEW (Fortran)" || EXIT_CODE=1
echo ""
compare_exact "$OUT_NEW" "$OUT_C" "new" "c" "NEW vs C" || EXIT_CODE=1
echo ""

# Statistical comparison
echo ""
echo "--- Statistical Verification ---"
if [ -f "compare_chains.R" ]; then
    echo "Checking REF vs NEW:"
    Rscript compare_chains.R "$OUT_REF/ref_train.hyp" "$OUT_NEW/new_train.hyp"
    [ $? -ne 0 ] && EXIT_CODE=1

    echo "Checking NEW vs C:"
    Rscript compare_chains.R "$OUT_NEW/new_train.hyp" "$OUT_C/c_train.hyp"
    [ $? -ne 0 ] && EXIT_CODE=1
else
    echo "Warning: compare_chains.R not found, skipping statistical check."
fi

if [ $EXIT_CODE -eq 0 ]; then
    echo "ALL VERIFICATIONS SUCCESSFUL"
else
    echo "VERIFICATION FAILED"
fi

exit $EXIT_CODE
