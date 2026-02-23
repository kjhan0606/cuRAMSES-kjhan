#!/bin/bash
# Ghost Zone Exchange Benchmark - Large grids (1024^3, 2048^3)
# Usage: sbatch -p a10 run_ghost_large.sh

#SBATCH -J ghost_large
#SBATCH -N 3
#SBATCH --ntasks-per-node=64
#SBATCH --exclude=syn01
#SBATCH -t 08:00:00
#SBATCH -o ghost_large_%j.log

EXE=/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/dd_speed/ghost_bench

echo "================================================"
echo " Ghost Zone Exchange Benchmark - Large Grids"
echo " Nodes: $SLURM_NNODES, Total slots: $SLURM_NTASKS"
echo " Date: $(date)"
echo "================================================"
echo ""

# N=1024, nvar=5, nrepeat=3
NGRID=1024
NVAR=5
NREPEAT=3
echo "########################################"
echo "# N=$NGRID  nrepeat=$NREPEAT"
echo "########################################"
for NP in 4 8 12 16 32 64; do
    echo "===== ncpu=$NP ====="
    mpirun -np $NP $EXE $NGRID $NVAR $NREPEAT
    echo ""
done

# N=2048, nvar=5, nrepeat=1
NGRID=2048
NVAR=5
NREPEAT=1
echo "########################################"
echo "# N=$NGRID  nrepeat=$NREPEAT"
echo "########################################"
for NP in 4 8 12; do
    echo "===== ncpu=$NP ====="
    mpirun -np $NP $EXE $NGRID $NVAR $NREPEAT
    echo ""
done

echo "Done: $(date)"
