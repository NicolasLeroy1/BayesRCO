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

# Output directories
OUT_REF="out_ref"
OUT_NEW="out_new"

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
        "../$bin" -bfile toy_app -out "${prefix}_train" -seed $SEED -ncat $NCAT -catfile toy_annot.txt -burnin $BURNIN -numit $NUMIT -thin $THIN
        "../$bin" -bfile toy_test -predict -out "${prefix}_test" -model "${prefix}_train.model" -freq "${prefix}_train.frq" -param "${prefix}_train.param" -ncat $NCAT -catfile toy_annot.txt
    )
}

compare_exact() {
    local dir1=$1
    local dir2=$2
    local prefix1=$3
    local prefix2=$4
    local diff_found=0

    echo "--- Comparing Exact Parity ($dir1 vs $dir2) ---"
    for ext in model param frq gv; do
        file1="$dir1/${prefix1}_train.$ext"
        [ "$ext" == "gv" ] && file1="$dir1/${prefix1}_test.gv"
        
        file2="$dir2/${prefix2}_train.$ext"
        [ "$ext" == "gv" ] && file2="$dir2/${prefix2}_test.gv"

        if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
            echo "Error: Missing files for .$ext comparison"
            diff_found=1
            continue
        fi

        if diff -q "$file1" "$file2" > /dev/null; then
            echo "[PASS] .$ext files are identical"
        else
            echo "[FAIL] .$ext files differ"
            diff_found=1
        fi
    done
    return $diff_found
}

# Main Execution
echo "Starting BayesRCO Verification System"

# Build
echo "Building versions..."
make clean
make ref new || { echo "Build failed"; exit 1; }

# Run versions
run_model "$BIN_REF" "$OUT_REF" "Reference (Old)" "ref"
run_model "$BIN_NEW" "$OUT_NEW" "New Modular (Fortran)" "new"

# Compare
EXIT_CODE=0

# Exact comparison between Fortran versions
compare_exact "$OUT_REF" "$OUT_NEW" "ref" "new" || EXIT_CODE=1

# Statistical comparison
echo "--- Statistical Verification ---"
if [ -f "compare_chains.R" ]; then
    Rscript compare_chains.R "$OUT_REF/ref_train.hyp" "$OUT_NEW/new_train.hyp"
    [ $? -ne 0 ] && EXIT_CODE=1
else
    echo "Warning: compare_chains.R not found, skipping statistical check."
fi

if [ $EXIT_CODE -eq 0 ]; then
    echo "VERIFICATION SUCCESSFUL"
else
    echo "VERIFICATION FAILED"
fi

exit $EXIT_CODE
