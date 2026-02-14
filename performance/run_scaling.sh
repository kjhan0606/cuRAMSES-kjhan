#!/bin/bash
# Scaling test: restart from HDF5 checkpoint and run 5 coarse steps
# Checkpoint: output_00002 (step 5, written with 12 CPUs)
# Variable-ncpu restart allows any CPU count
# Usage: bash run_scaling.sh

BASEDIR="$(cd "$(dirname "$0")" && pwd)"
EXE="$BASEDIR/../bin/ramses_final3d"
CHECKPOINT="$BASEDIR/../test_ksection/run_varcpu_12cpu/output_00002"
YIELD="$BASEDIR/../test_ksection/yield_table.asc"

# CPU counts to test
CPUS="4 8 12 16 24 32 48 64"

# Check prerequisites
if [ ! -f "$EXE" ]; then
    echo "ERROR: Executable not found: $EXE"
    exit 1
fi
if [ ! -d "$CHECKPOINT" ]; then
    echo "ERROR: Checkpoint not found: $CHECKPOINT"
    exit 1
fi

cd "$BASEDIR"

for N in $CPUS; do
    RUNDIR="run_${N}cpu"
    echo "============================================"
    echo "  Running with $N CPUs (restart from checkpoint)"
    echo "============================================"

    # Compute ngridtot/nparttot: ensure ngridmax >= 3M per rank
    NGRIDTOT=$((N * 3000000))
    if [ $NGRIDTOT -lt 40000000 ]; then
        NGRIDTOT=40000000
    fi
    NPARTTOT=$((N * 4000000))
    if [ $NPARTTOT -lt 200000000 ]; then
        NPARTTOT=200000000
    fi

    mkdir -p $RUNDIR
    cd $RUNDIR

    # Generate per-run namelist
    cat > cosmo_perf.nml << ENDNML
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

    # Link yield table and checkpoint
    ln -sf "$YIELD" yield_table.asc
    ln -sf "$CHECKPOINT" output_00002

    echo "  ngridtot=$NGRIDTOT nparttot=$NPARTTOT"

    mpirun -np $N "$EXE" cosmo_perf.nml > perf_${N}cpu.log 2>&1
    EXIT=$?

    if [ $EXIT -ne 0 ]; then
        echo "  FAILED (exit code $EXIT)"
    else
        # Extract timing info (Total elapsed time is on a single line)
        TOTAL=$(grep "Total elapsed time:" perf_${N}cpu.log | awk '{print $4}')
        # Main step format: Main step=     10 mcons= X econs= Y epot=Z ekin= W
        LINE=$(grep "Main step=     10" perf_${N}cpu.log)
        ECONS=$(echo "$LINE" | sed 's/.*econs= *//' | awk '{print $1}')
        EPOT=$(echo "$LINE" | sed 's/.*epot=//' | awk '{print $1}')
        EKIN=$(echo "$LINE" | sed 's/.*ekin= *//' | awk '{print $1}')
        echo "  Elapsed: ${TOTAL}s"
        echo "  econs=$ECONS epot=$EPOT ekin=$EKIN"
    fi

    cd "$BASEDIR"
    echo ""
done

# Summary table
echo ""
echo "============================================================"
echo "  SCALING TEST SUMMARY"
echo "  Checkpoint restart (step 5 -> 10, 5 coarse steps)"
echo "============================================================"
printf "%-6s %-14s %-12s %-14s %-14s\n" "ncpu" "elapsed(s)" "econs" "epot" "ekin"
echo "----------------------------------------------------------------------"
for N in $CPUS; do
    LOGFILE="run_${N}cpu/perf_${N}cpu.log"
    if [ -f "$LOGFILE" ]; then
        TOTAL=$(grep "Total elapsed time:" "$LOGFILE" | awk '{print $4}')
        LINE=$(grep "Main step=     10" "$LOGFILE")
        if [ -z "$LINE" ]; then
            printf "%-6s %-14s\n" "$N" "FAILED"
        else
            ECONS=$(echo "$LINE" | sed 's/.*econs= *//' | awk '{print $1}')
            EPOT=$(echo "$LINE" | sed 's/.*epot=//' | awk '{print $1}')
            EKIN=$(echo "$LINE" | sed 's/.*ekin= *//' | awk '{print $1}')
            printf "%-6s %-14s %-12s %-14s %-14s\n" "$N" "$TOTAL" "$ECONS" "$EPOT" "$EKIN"
        fi
    else
        printf "%-6s %-14s\n" "$N" "NO LOG"
    fi
done

# Detailed timer breakdown
echo ""
echo "============================================================"
echo "  TIMER BREAKDOWN (average seconds per rank)"
echo "============================================================"
printf "%-6s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n" \
    "ncpu" "coarse" "refine" "loadbal" "particle" "poisson" "poi-mg" "hydro-gz" "godunov"
echo "----------------------------------------------------------------------------------------------"
for N in $CPUS; do
    LOGFILE="run_${N}cpu/perf_${N}cpu.log"
    if [ -f "$LOGFILE" ] && grep -q "Main step=     10" "$LOGFILE"; then
        COARSE=$(grep "coarse levels" "$LOGFILE" | awk '{print $2}')
        REFINE=$(grep "refine " "$LOGFILE" | head -1 | awk '{print $2}')
        LOADBAL=$(grep "loadbalance" "$LOGFILE" | awk '{print $2}')
        PARTICLE=$(grep "particles " "$LOGFILE" | head -1 | awk '{print $2}')
        POISSON=$(grep "poisson " "$LOGFILE" | head -1 | awk '{print $2}')
        POISMG=$(grep "poisson - mg" "$LOGFILE" | awk '{print $2}')
        HYDROGZ=$(grep "hydro - ghostzones" "$LOGFILE" | head -1 | awk '{print $2}')
        GODUNOV=$(grep "hydro - godunov" "$LOGFILE" | awk '{print $2}')
        printf "%-6s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n" \
            "$N" "$COARSE" "$REFINE" "$LOADBAL" "$PARTICLE" "$POISSON" "$POISMG" "$HYDROGZ" "$GODUNOV"
    fi
done
