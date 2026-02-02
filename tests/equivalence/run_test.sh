#!/bin/bash
# Multi-chain equivalence verification test for BayesRCO implementations
# Runs 5 chains (seeds 41-45) for REF, NEW, and C implementations
# Compares distributions across chains to verify statistical equivalence

set -e

# Navigate to project root
cd "$(dirname "$0")/../.."
ROOT_DIR=$(pwd)

# Configuration
BIN_DIR="bin"
TOY_DIR="toy_dataset"
NUMIT=300
THIN=5
BURNIN=0
NCAT=4
SEEDS=(41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60)

# Binaries
BIN_REF="$BIN_DIR/bayesRCO_ref"
BIN_NEW="$BIN_DIR/bayesRCO"
BIN_C="$BIN_DIR/bayesRCO_c"

# Output directories
OUT_DIR="tests/equivalence/outputs"

# Functions
run_chain() {
    local bin=$1
    local out_dir=$2
    local seed=$3
    
    mkdir -p "$out_dir"
    cp $TOY_DIR/toy_annot.txt $TOY_DIR/toy_app.* $TOY_DIR/toy_test.* "$out_dir/"
    
    (
        cd "$out_dir" || exit 1
        "$ROOT_DIR/$bin" -bfile toy_app -out train -seed $seed -ncat $NCAT -catfile toy_annot.txt -burnin $BURNIN -numit $NUMIT -thin $THIN > /dev/null 2>&1
    )
}

# Main Execution
echo "=== BayesRCO Multi-Chain Equivalence Verification ==="
echo "Seeds: ${SEEDS[*]}"
echo "Working directory: $ROOT_DIR"
echo ""

# Build using centralized Makefile
echo "Building all versions..."
make clean > /dev/null 2>&1
make all > /dev/null 2>&1 || { echo "Build failed"; exit 1; }
echo "Build complete."
echo ""

# Run chains
echo "Running chains..."
for seed in "${SEEDS[@]}"; do
    echo -n "  Seed $seed: REF..."
    run_chain "$BIN_REF" "$OUT_DIR/ref_$seed" $seed
    echo -n " NEW..."
    run_chain "$BIN_NEW" "$OUT_DIR/new_$seed" $seed
    echo -n " C..."
    run_chain "$BIN_C" "$OUT_DIR/c_$seed" $seed
    echo " done"
done
echo ""

# Combine hyp files for comparison
echo "Combining chain outputs..."
rm -f "$OUT_DIR/ref_combined.hyp" "$OUT_DIR/new_combined.hyp" "$OUT_DIR/c_combined.hyp"
for seed in "${SEEDS[@]}"; do
    # Skip header lines, add data with seed column
    awk -v s="$seed" 'NR>1 {print $0, s}' "$OUT_DIR/ref_$seed/train.hyp" >> "$OUT_DIR/ref_combined.hyp" 2>/dev/null || true
    awk -v s="$seed" 'NR>1 {print $0, s}' "$OUT_DIR/new_$seed/train.hyp" >> "$OUT_DIR/new_combined.hyp" 2>/dev/null || true
    awk -v s="$seed" 'NR>1 {print $0, s}' "$OUT_DIR/c_$seed/train.hyp" >> "$OUT_DIR/c_combined.hyp" 2>/dev/null || true
done
echo "Combined files created."
echo ""

# Statistical comparison
echo "--- Statistical Verification (5 chains each) ---"
COMPARE_SCRIPT="tests/equivalence/compare_chains.R"

echo ""
echo "Checking REF vs NEW:"
Rscript "$COMPARE_SCRIPT" "$OUT_DIR/ref_combined.hyp" "$OUT_DIR/new_combined.hyp" "$OUT_DIR/ref_vs_new.pdf"
REF_NEW_STATUS=$?

echo ""
echo "Checking NEW vs C:"
Rscript "$COMPARE_SCRIPT" "$OUT_DIR/new_combined.hyp" "$OUT_DIR/c_combined.hyp" "$OUT_DIR/new_vs_c.pdf"
NEW_C_STATUS=$?

echo ""
if [ $REF_NEW_STATUS -eq 0 ] && [ $NEW_C_STATUS -eq 0 ]; then
    echo "=== ALL VERIFICATIONS SUCCESSFUL ==="
    exit 0
else
    echo "=== VERIFICATION COMPLETED WITH DIFFERENCES ==="
    exit 1
fi
