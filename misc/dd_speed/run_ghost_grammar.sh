#!/bin/bash
# Ghost Zone Exchange Benchmark - grammar cluster max config
# 16 nodes x 64 CPUs = 1024 ranks, N=1024 and N=2048

#SBATCH -J ghost_max
#SBATCH -p normal
#SBATCH -N 16
#SBATCH --ntasks-per-node=64
#SBATCH -t 04:00:00
#SBATCH -o ghost_max_%j.log

EXE=/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/dd_speed/ghost_bench

echo "================================================"
echo " Ghost Zone Exchange Benchmark - grammar MAX"
echo " Nodes: $SLURM_NNODES, Total ranks: $SLURM_NTASKS"
echo " Date: $(date)"
echo "================================================"
echo ""

# --- N=1024, scaling from 64 to 1024 ranks ---
NGRID=1024
NVAR=5
NREPEAT=3
echo "########################################"
echo "# N=$NGRID  nrepeat=$NREPEAT"
echo "########################################"
for NP in 64 128 256 512 1024; do
    echo ""
    echo "===== ncpu=$NP ====="
    mpirun -np $NP $EXE $NGRID $NVAR $NREPEAT
done

echo ""
echo ""

# --- N=2048, max ranks ---
NGRID=2048
NVAR=5
NREPEAT=1
echo "########################################"
echo "# N=$NGRID  nrepeat=$NREPEAT"
echo "########################################"
for NP in 256 512 1024; do
    echo ""
    echo "===== ncpu=$NP ====="
    mpirun -np $NP $EXE $NGRID $NVAR $NREPEAT
done

echo ""
echo "Done: $(date)"
