#!/bin/bash
# Strong scaling test: MPI ranks × OpenMP threads
# Checkpoint: output_00002 (step 5 → 10, 5 coarse steps)
# Machine: 2× AMD EPYC 7543 (64 physical cores, 128 logical)
# Usage: bash run_strong_scaling.sh

set -e
BASEDIR="$(cd "$(dirname "$0")" && pwd)"
EXE="$BASEDIR/../bin/ramses_final3d"
CHECKPOINT="$BASEDIR/../test_ksection/run_varcpu_12cpu/output_00002"
YIELD="$BASEDIR/../test_ksection/yield_table.asc"
OUTDIR="$BASEDIR/strong_scaling"
mkdir -p "$OUTDIR"

# Test configurations: "nranks:nthreads" pairs
# Group 1: Pure MPI (1 thread), varying ranks
# Group 2: Fixed 64 cores, varying MPI/OMP ratio
CONFIGS=(
    "2:1"
    "4:1"
    "8:1"
    "16:1"
    "32:1"
    "64:1"
    "32:2"
    "16:4"
    "8:8"
    "4:16"
)

run_one() {
    local NR=$1   # MPI ranks
    local NT=$2   # OMP threads
    local TAG="r${NR}_t${NT}"
    local RUNDIR="$OUTDIR/$TAG"
    local LOGFILE="$RUNDIR/run.log"

    echo "--- $TAG (${NR} ranks × ${NT} threads = $((NR*NT)) cores) ---"

    # Compute ngridtot/nparttot per rank
    local NGRIDTOT=$((NR * 3000000))
    [ $NGRIDTOT -lt 40000000 ] && NGRIDTOT=40000000
    local NPARTTOT=$((NR * 20000000))
    [ $NPARTTOT -lt 200000000 ] && NPARTTOT=200000000

    mkdir -p "$RUNDIR"

    cat > "$RUNDIR/cosmo.nml" << ENDNML
&RUN_PARAMS
cosmo=.true.
pic=.true.
poisson=.true.
hydro=.true.
nrestart=2
nremap=5
nsubcycle=1,1,2
ncontrol=1
nstepmax=10
ordering='ksection'
memory_balance=.true.
/

&OUTPUT_PARAMS
noutput=1
aout=1.0
foutput=10
informat='hdf5'
outformat='hdf5'
/

&INIT_PARAMS
filetype='grafic'
initfile(1)='/gpfs/kjhan/Darwin/Cudaization/MUSIC/ics_ramses/level_008'
/

&AMR_PARAMS
levelmin=8
levelmax=10
nexpand=1
ngridtot=${NGRIDTOT}
nparttot=${NPARTTOT}
/

&LIGHTCONE_PARAMS
/

&REFINE_PARAMS
m_refine=3*8.,
ivar_refine=0
interpol_var=1
interpol_type=0
/

&HYDRO_PARAMS
gamma=1.6666667
courant_factor=0.8
scheme='muscl'
slope_type=1
/

&POISSON_PARAMS
/

&PHYSICS_PARAMS
sf_birth_properties=.false.
yieldtablefilename='yield_table.asc'
/
ENDNML

    ln -sf "$YIELD" "$RUNDIR/yield_table.asc"
    ln -sf "$CHECKPOINT" "$RUNDIR/output_00002"

    # Set OpenMP environment
    export OMP_NUM_THREADS=$NT
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores

    # Run from the run directory (HDF5 uses relative paths)
    pushd "$RUNDIR" > /dev/null
    mpirun -np $NR -genv I_MPI_PIN_DOMAIN=omp --bind-to none \
        "$EXE" cosmo.nml > run.log 2>&1
    local RC=$?
    popd > /dev/null

    if [ $RC -ne 0 ]; then
        echo "  FAILED (exit $RC)"
        return 1
    fi

    # Extract results
    local TOTAL=$(grep "Total elapsed time:" "$LOGFILE" | awk '{print $4}')
    local LINE=$(grep "Main step=     10" "$LOGFILE")
    local ECONS=$(echo "$LINE" | sed 's/.*econs= *//' | awk '{print $1}')
    local EPOT=$(echo "$LINE" | sed 's/.*epot=//' | awk '{print $1}')
    local EKIN=$(echo "$LINE" | sed 's/.*ekin= *//' | awk '{print $1}')
    echo "  Time: ${TOTAL}s  econs=$ECONS epot=$EPOT ekin=$EKIN"
}

echo "========================================================"
echo "  STRONG SCALING TEST"
echo "  Problem: 200M particles, levelmin=8, levelmax=10"
echo "  Machine: 2× EPYC 7543 (64 cores), 1 TB RAM"
echo "  Steps: 5 coarse steps (restart step 5→10)"
echo "========================================================"
echo ""

for CFG in "${CONFIGS[@]}"; do
    NR=${CFG%%:*}
    NT=${CFG##*:}
    run_one $NR $NT
    echo ""
done

# Generate summary table
echo ""
echo "========================================================"
echo "  SUMMARY"
echo "========================================================"
printf "%-8s %-8s %-8s %-12s %-12s %-14s %-14s\n" \
    "nranks" "nthrd" "ncores" "elapsed(s)" "speedup" "epot" "ekin"
echo "------------------------------------------------------------------------"

# Get baseline (2 ranks × 1 thread) time for speedup calculation
BASELINE=""
for CFG in "${CONFIGS[@]}"; do
    NR=${CFG%%:*}
    NT=${CFG##*:}
    TAG="r${NR}_t${NT}"
    LOGFILE="$OUTDIR/$TAG/run.log"
    if [ -f "$LOGFILE" ] && grep -q "Main step=     10" "$LOGFILE"; then
        TOTAL=$(grep "Total elapsed time:" "$LOGFILE" | awk '{print $4}')
        if [ -z "$BASELINE" ]; then
            BASELINE=$TOTAL
        fi
        SPEEDUP=$(echo "scale=2; $BASELINE / $TOTAL" | bc 2>/dev/null || echo "N/A")
        LINE=$(grep "Main step=     10" "$LOGFILE")
        EPOT=$(echo "$LINE" | sed 's/.*epot=//' | awk '{print $1}')
        EKIN=$(echo "$LINE" | sed 's/.*ekin= *//' | awk '{print $1}')
        printf "%-8s %-8s %-8s %-12s %-12s %-14s %-14s\n" \
            "$NR" "$NT" "$((NR*NT))" "$TOTAL" "${SPEEDUP}x" "$EPOT" "$EKIN"
    else
        printf "%-8s %-8s %-8s %-12s\n" "$NR" "$NT" "$((NR*NT))" "FAILED"
    fi
done

# Timer breakdown
echo ""
echo "========================================================"
echo "  TIMER BREAKDOWN (average seconds)"
echo "========================================================"
printf "%-8s %-8s %-9s %-9s %-9s %-9s %-9s %-9s %-9s %-9s\n" \
    "nranks" "nthrd" "coarse" "refine" "loadbal" "particle" "poisson" "poi-mg" "hydro-gz" "godunov"
echo "----------------------------------------------------------------------------------------------"
for CFG in "${CONFIGS[@]}"; do
    NR=${CFG%%:*}
    NT=${CFG##*:}
    TAG="r${NR}_t${NT}"
    LOGFILE="$OUTDIR/$TAG/run.log"
    if [ -f "$LOGFILE" ] && grep -q "Main step=     10" "$LOGFILE"; then
        COARSE=$(grep "coarse levels" "$LOGFILE" | awk '{print $2}')
        REFINE=$(grep "refine " "$LOGFILE" | head -1 | awk '{print $2}')
        LOADBAL=$(grep "loadbalance" "$LOGFILE" | awk '{print $2}')
        PARTICLE=$(grep "particles " "$LOGFILE" | head -1 | awk '{print $2}')
        POISSON=$(grep "^  *[0-9].*poisson " "$LOGFILE" | head -1 | awk '{print $2}')
        POISMG=$(grep "poisson - mg" "$LOGFILE" | awk '{print $2}')
        HYDROGZ=$(grep "hydro - ghostzones" "$LOGFILE" | head -1 | awk '{print $2}')
        GODUNOV=$(grep "hydro - godunov" "$LOGFILE" | awk '{print $2}')
        printf "%-8s %-8s %-9s %-9s %-9s %-9s %-9s %-9s %-9s %-9s\n" \
            "$NR" "$NT" "$COARSE" "$REFINE" "$LOADBAL" "$PARTICLE" "$POISSON" "$POISMG" "$HYDROGZ" "$GODUNOV"
    fi
done
