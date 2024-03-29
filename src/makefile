
# MESRA makefile

.SUFFIXES: .o .f90

# Selects the compiler.

#FC=ifort  
FC=gfortran

# Chooses the compilation parameters.

### For development and bug tracking

# Used when testing gfortran
#FFLAGS=-g -I. -O0 -c 

# Used when testing ifort
#FFLAGS = -O0 -C -g -c

### For execution

# With gfortran
FFLAGS=-g -I. -O2 -c 

# With ifort
#FFLAGS = -O2 -g -c

# Sets the libraries and their path.

# intel libraries
#LFLAGS= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread 
#LFLAGS = -mkl=parallel

# blas and lapack libraries
#LFLAGS= -llapack -lblas
#LFLAGS=-L/opt/lapack/gnu -llapack -lblas
LFLAGS=-L/usr/lib -llapack -lblas

# Defines a .o for any file entering the compilation.

OBJS= \
declare.o rgen_info.o rPA_status.o rdaggerf.o rinp.o \
rfchk.o orth.o rC.o rS.o rP.o rP_g.o rPx_g.o rPxspin_g.o rP_qp.o \
mo_to_ao.o ao_to_mo.o double_mat_prod.o rPx_qp.o \
trace_mat.o det_at.o weigen.o dau_main.o \
zvec_main.o daz_main.o da.o launch_svd_and_rotate.o svd_and_rotation_XY.o svd.o diag.o \
print_mat_mo_to_ao_fchk.o \
triangleK.o PopulationAnalysis.o  wfchk_dens.o wfchk_orbsAlpha.o wfchk_orbsBeta.o \
PA_mat.o QuantumDescriptorsPA.o rPrelax_g.o PA_analysis.o \
common_dar_daz.o adiab_z.o no_adiab_z.o \
rTK.o common_XY.o daXY.o da_XY.o trans_orbs_XY.o orbs_XY.o \
deal.o increase_cube_size.o \
rlxy_PA.o rlxy_PAops.o \
rqm_NI_info.o rdcubehdr.o consistency.o rddensities.o printingdadens.o qm_NI.o \
split.o cubeop.o \
alphaddagger.o alpha_ddagger.o qmNIrlx.o \
main_mesra.o \

# Command line for the creation of the .o files.

.f90.o:

	$(FC) $(FFLAGS) $*.f90

# Creation of the main file.

all: main_mesra

main_mesra :     $(OBJS)
	$(FC) $(OBJS)  $(LFLAGS) -o main_mesra

# Removal of the .o files after their compilation into the main file.

clean : 
	rm -f *.o *.mod

