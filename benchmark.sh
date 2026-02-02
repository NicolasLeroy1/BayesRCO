#!/bin/bash

# Configuration
TOY_DIR="toy_dataset"
OUT_BM="outputs/benchmark"
NUMIT=5000
BURNIN=0

BIN_REF="bin/bayesRCO_ref"
BIN_NEW="bin/bayesRCO"
BIN_C="bin/bayesRCO_c"

mkdir -p $OUT_BM
cp $TOY_DIR/toy_annot.txt $TOY_DIR/toy_app.* $TOY_DIR/toy_test.* $OUT_BM/

echo "----------------------------------------------------------"
echo "Benchmarking BayesRCO Versions (numit=$NUMIT, burnin=$BURNIN)"
echo "----------------------------------------------------------"

# Helper function to run and time
run_benchmark() {
    NAME=$1
    CMD=$2
    LOG="${NAME}.log"
    TIME_LOG="${NAME}_time.txt"
    
    echo "Running $NAME..."
    /usr/bin/time -v $CMD > $LOG 2> $TIME_LOG
    
    # Extract stats
    USER_TIME=$(grep "User time" $TIME_LOG | cut -d: -f2 | xargs)
    SYS_TIME=$(grep "System time" $TIME_LOG | cut -d: -f2 | xargs)
    WALL_TIME=$(grep "Elapsed (wall clock)" $TIME_LOG | cut -d: -f2- | xargs)
    MAX_MEM=$(grep "Maximum resident set size" $TIME_LOG | cut -d: -f2 | xargs)
    
    echo "  Time (User): $USER_TIME s"
    echo "  Time (Sys):  $SYS_TIME s"
    echo "  Wall Time:   $WALL_TIME"
    echo "  Max Memory:  $MAX_MEM kbytes"
    echo ""
}

# Go to benchmark directory
# We use absolute paths for binaries to avoid ../ issues if moved
ROOT_DIR=$(pwd)
cd $OUT_BM

# 1. Reference Fortran
run_benchmark "Reference_Fortran" "$ROOT_DIR/$BIN_REF -bfile toy_app -out toy_bm_ref -seed 10 -ncat 4 -catfile toy_annot.txt -burnin $BURNIN -numit $NUMIT"

# 2. New Fortran
run_benchmark "New_Fortran" "$ROOT_DIR/$BIN_NEW -bfile toy_app -out toy_bm_new -seed 10 -ncat 4 -catfile toy_annot.txt -burnin $BURNIN -numit $NUMIT"

# 3. C Version
run_benchmark "C_Version" "$ROOT_DIR/$BIN_C -bfile toy_app -out toy_bm_c -seed 10 -ncat 4 -catfile toy_annot.txt -burnin $BURNIN -numit $NUMIT"

cd $ROOT_DIR
