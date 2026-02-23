FC = mpiifx
FFLAGS = -O2 -xHost

all: dd_bench ghost_bench mpi_bench lb_bench poisson_bench

dd_bench: dd_bench.f90
	$(FC) $(FFLAGS) -o $@ $<

ghost_bench: ghost_bench.f90
	$(FC) $(FFLAGS) -o $@ $<

mpi_bench: mpi_bench.f90
	$(FC) $(FFLAGS) -o $@ $<

lb_bench: lb_bench.f90
	$(FC) $(FFLAGS) -o $@ $<

poisson_bench: poisson_bench.f90
	$(FC) $(FFLAGS) -o $@ $<

clean:
	rm -f dd_bench ghost_bench mpi_bench lb_bench poisson_bench *.o *.mod
