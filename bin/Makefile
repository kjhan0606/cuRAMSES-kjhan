#############################################################################
# If you have problems with this makefile, contact Romain.Teyssier@gmail.com
#############################################################################
# Compilation time parameters
#NVECTOR = 32
#NDIM = 3
#NPRE = 8
#NVAR = 16
#NENER = 0
#SOLVER = hydro
#PATCH = ../patch/Horizon5-master
#EXEC = ramses_ch
##################################
# Yohan's patch
#NVECTOR = 32
NVECTOR = 32
NDIM = 3
NPRE = 8
NVAR = 11
NENER = 0
SOLVER = hydro
PATCH = ../patch/Horizon5-master-2
EXEC = ramses_final
#############################################################################
GITBRANCH = $(shell git rev-parse --abbrev-ref HEAD)
GITHASH = $(shell git log --pretty=format:'%H' -n 1)
GITREPO = $(shell git config --get remote.origin.url)
BUILDDATE = $(shell date +"%D-%T")
DEFINES = -DNVECTOR=$(NVECTOR) -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DNENER=$(NENER) -DNVAR=$(NVAR) \
          -DSOLVER=$(SOLVER)
#############################################################################
# Fortran compiler options and directives

# --- MPI, pgf90 syntax ------------------------------
#F90 = mpif90 -O3
#FFLAGS = -Mpreprocess $(DEFINES) -Mbackslash

# --- MPI, ifort syntax ------------------------------
#F90 = mpif90 -openmp #-xMIC-AVX512
F90 = mpiifx -qopenmp -march=core-avx2
FFLAGS = -cpp -O3 $(DEFINES) -DNOSYSTEM -DLONGINT -DQUADHILBERT -DOUTPUT_PARTICLE_POTENTIAL

#############################################################################
# HDF5 parallel I/O support (build with: make HDF5=1)
#############################################################################
ifdef HDF5
  HDF5_DIR = $(HOME)/local/hdf5
  FFLAGS += -DHDF5 -I$(HDF5_DIR)/include
  LIBS_HDF5 = -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 -lz -Wl,-rpath,$(HDF5_DIR)/lib
  HDF5_MODOBJ = ramses_hdf5_io.o
  HDF5_OBJS = backup_hdf5.o restore_hdf5.o
else
  LIBS_HDF5 =
  HDF5_MODOBJ =
  HDF5_OBJS =
endif

#############################################################################
# CUDA GPU hydro support (build with: make USE_CUDA=1)
#############################################################################
ifdef USE_CUDA
  CUDA_DIR = /opt/ohpc/pub/cuda/13.0.2
  NVCC = $(CUDA_DIR)/bin/nvcc
  CUDA_ARCH = -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_89,code=sm_89 -gencode arch=compute_90,code=sm_90
  N_CUDA_STREAMS ?= 1
  NVCCFLAGS = -O3 $(CUDA_ARCH) -DNDIM=$(NDIM) -DNVAR=$(NVAR) -DNVECTOR=$(NVECTOR) -DN_STREAMS=$(N_CUDA_STREAMS)
  FFLAGS += -DHYDRO_CUDA
  LIBS_CUDA = -L$(CUDA_DIR)/lib64 -lcudart -lcufft -lstdc++ -Wl,-rpath,$(CUDA_DIR)/lib64
  CUDA_MODOBJ = cuda_commons.o hydro_cuda_interface.o poisson_cuda_interface.o
  CUDA_OBJS = cuda_stream_pool.o hydro_cuda_kernels.o poisson_cuda_kernels.o
else
  LIBS_CUDA =
  CUDA_MODOBJ =
  CUDA_OBJS =
endif

#############################################################################
MOD = mod
#############################################################################
# MPI librairies
LIBMPI =
#LIBMPI = -lfmpi -lmpi -lelan

# --- CUDA libraries, for Titane ---
LIBCUDA = -L/opt/cuda/lib  -lm -lcuda -lcudart

LIBGRACKLE = -lgrackle -lhdf5 -lz
LIBS = $(LIBMPI) $(LIBS_HDF5) $(LIBS_CUDA)
#############################################################################
# Sources directories are searched in this exact order
VPATH = ../patch/cuda:../patch/oct_tree:$(PATCH):../$(SOLVER):../aton:../hydro:../pm:../poisson:../amr
#############################################################################
# All objects
MODOBJ	= amr_parameters.jaehyun.o amr_commons.kjhan.o random.o pm_parameters.o pm_commons.o poisson_parameters.o \
	 poisson_commons.o hydro_parameters.o hydro_commons.o cooling_module.o bisection.o ksection.o sparse_mat.o \
	 morton_keys.o morton_hash.o \
	 clfind_commons.o gadgetreadfile.o write_makefile.o write_gitinfo.o \
	 $(HDF5_MODOBJ) $(CUDA_MODOBJ)

AMROBJ = read_params.jaehyun.o init_amr.o init_time.o init_refine.o adaptive_loop.jaehyun.o amr_step.jaehyun.o update_time.o \
         output_amr.kjhan.o flag_utils.kjhan.o physical_boundaries.kjhan.o virtual_boundaries.kjhan.o refine_utils.o nbors_utils.kjhan.o \
         hilbert.o load_balance.kjhan.o title.o sort.o cooling_fine.kjhan.o units.o light_cone.part.o light_cone.hydro2.o light_cone.sink2.o movie_mod.yonghwi_org.o kjhan.o\
		 output_sphere_hydro.o output_sphere_part.o output_sphere_sink.o \
         morton_init.o $(HDF5_OBJS)
# Particle-Mesh objects
PMOBJ = init_part.o output_part.o rho_fine.kjhan.o synchro_fine.kjhan.o move_fine.o newdt_fine.kjhan.o particle_tree.kjhan.o \
        add_list.o remove_list.o star_formation.kjhan.o sink_particle.kjhan.o feedback.kjhan3.o \
		clump_finder.o clump_merger.o \
        flag_formation_sites.o init_sink.o output_sink.o
# Poisson solver objects
#POISSONOBJ = init_poisson.o phi_fine_cg.kisti.o interpol_phi.kjhan.o force_fine.kjhan.o multigrid_coarse.kjhan.o multigrid_fine_commons.o \
             multigrid_fine_fine.kisti.o multigrid_fine_coarse.kisti.o gravana.o boundary_potential.o rho_ana.kjhan.o output_poisson.o
POISSONOBJ = init_poisson.o phi_fine_cg.kjhan.o interpol_phi.kjhan.o force_fine.kjhan.o multigrid_coarse.kjhan.o multigrid_fine_commons.o \
             multigrid_fine_fine.kjhan.o multigrid_fine_coarse.kjhan.o gravana.o boundary_potential.kjhan.o rho_ana.kjhan.o output_poisson.o
# Hydro objects
HYDROOBJ = init_hydro.o init_flow_fine.o write_screen.o output_hydro.o courant_fine.kjhan.o godunov_fine.kjhan.o \
           uplmde.kjhan.o umuscl.kjhan.o interpol_hydro.kjhan.o godunov_utils.kjhan.o condinit.o hydro_flag.kjhan.o hydro_boundary.o \
           boundana.o read_hydro_params.o synchro_hydro_fine.kjhan.o
# All objects
AMRLIB = $(AMROBJ) $(HYDROOBJ) $(PMOBJ) $(POISSONOBJ)
# ATON objects
ATON_MODOBJ = timing.o radiation_commons.o rad_step.o
ATON_OBJ = observe.o init_radiation.o rad_init.o rad_boundary.o rad_stars.o rad_backup.o ../aton/atonlib/libaton.a
#############################################################################
ramses:	$(MODOBJ) $(AMRLIB) $(CUDA_OBJS) ramses.o
	$(F90) $(MODOBJ) $(AMRLIB) $(CUDA_OBJS) ramses.o -o $(EXEC)$(NDIM)d $(LIBS)
	rm write_makefile.f90
#	rm write_patch.f90
ramses_aton: $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o
	$(F90) $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o -o $(EXEC)$(NDIM)d $(LIBS) $(LIBCUDA)
	rm write_makefile.f90
#	rm write_patch.f90
#############################################################################
write_gitinfo.o: FORCE
	$(F90) $(FFLAGS) -DPATCH=\'$(PATCH)\' -DGITBRANCH=\'$(GITBRANCH)\' -DGITHASH=\'"$(GITHASH)"\' \
 -DGITREPO=\'$(GITREPO)\' -DBUILDDATE=\'"$(BUILDDATE)"\' -c ../amr/write_gitinfo.f90 -o $@
write_makefile.o: FORCE
	../utils/scripts/cr_write_makefile.sh $(MAKEFILE_LIST)
	$(F90) $(FFLAGS) -c write_makefile.f90 -o $@
write_patch.o: FORCE
	../utils/scripts/cr_write_patch.sh $(PATCH)
	$(F90) $(FFLAGS) -c write_patch.f90 -o $@
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@
ifdef USE_CUDA
%.o:%.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
hydro_cuda_kernels.o: cuda_stream_pool.o
poisson_cuda_kernels.o: cuda_stream_pool.o
hydro_cuda_interface.o: cuda_commons.o
poisson_cuda_interface.o: cuda_commons.o
cuda_commons.o: amr_commons.kjhan.o
endif
FORCE:
#############################################################################
clean :
	rm -f *.o *.$(MOD)
#############################################################################
