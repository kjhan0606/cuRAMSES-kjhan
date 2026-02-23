#!/bin/bash
# Combined Benchmark Suite - grammar cluster max config
# 16 nodes x 64 CPUs = 1024 ranks

#SBATCH -J bench_all
#SBATCH -p normal
#SBATCH -N 16
#SBATCH --ntasks-per-node=64
#SBATCH -t 04:00:00
#SBATCH -o bench_all_%j.log

DIR=/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/dd_speed

echo "================================================"
echo " MPI Benchmark Suite - grammar MAX"
echo " Nodes: $SLURM_NNODES, Total ranks: $SLURM_NTASKS"
echo " Date: $(date)"
echo "================================================"
echo ""

#===========================================
# 1. MPI Communication Benchmark
#===========================================
echo "###################################################"
echo "# MPI Communication Benchmark"
echo "###################################################"
NREPEAT=100
for NP in 64 128 256 512 1024; do
    echo ""
    echo "===== mpi_bench ncpu=$NP nrepeat=$NREPEAT ====="
    mpirun -np $NP $DIR/mpi_bench $NREPEAT
done

echo ""
echo ""

#===========================================
# 2. Load Balance Benchmark
#===========================================
echo "###################################################"
echo "# Load Balance Benchmark"
echo "###################################################"
NCLUSTERS=20
for NGRID in 128 256 512; do
    for NP in 64 128 256 512 1024; do
        echo ""
        echo "===== lb_bench N=$NGRID ncpu=$NP nclusters=$NCLUSTERS ====="
        mpirun -np $NP $DIR/lb_bench $NGRID $NCLUSTERS
    done
done

echo ""
echo ""

#===========================================
# 3. Poisson Solver Benchmark
#===========================================
echo "###################################################"
echo "# Poisson Solver Benchmark"
echo "###################################################"
MAX_ITER=100
for NGRID in 128 256 512 1024; do
    # ncpu must divide N
    for NP in 64 128 256 512 1024; do
        # Check N mod ncpu == 0
        if [ $((NGRID % NP)) -eq 0 ]; then
            echo ""
            echo "===== poisson_bench N=$NGRID ncpu=$NP max_iter=$MAX_ITER ====="
            mpirun -np $NP $DIR/poisson_bench $NGRID $MAX_ITER
        fi
    done
done

echo ""
echo "Done: $(date)"
