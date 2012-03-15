###################################################################
#  Makefile for unigrid resistive MHD code
#  Daniel R. Reynolds, SMU Mathematics
###################################################################

# Automatic compilation of .o files from Fortran source
.SUFFIXES: .F .o

# Make.config is used for setting compile-time problem 
# configuration settings (boundary conditions, etc).
include Make.config

# Make.machine is used for setting machine-dependent compilation
# variables: e.g. compilers, libraries, F90-to-C conversion vars
include Make.machine

# append problem and interface preprocessor defs to compilation defs.
FFLAGS_ = $(FFLAGS) $(PROB_DEFS) $(SUND_DEFS) $(SILO_DEFS)
CFLAGS_ = $(CFLAGS) $(PROB_DEFS) $(SUND_DEFS) $(SILO_DEFS)



####### Executable source dependencies #######

# all executables depend on these files
ALLSRC = modules profiling main setupgrid setuplocaldomain \
         setboundaryvalues output mkdir_output write_silo  \
         bdryexchange seteigenvectors_cons inviscidfluxcd2 \
         inviscidfluxcd4 inviscidfluxtcd projection newdt  \
         xchange diags drivertools mhdtools derf           \
         seteigenvectors_prim constructlrstates4 setalphas \
         setdurl setsls setroevariables seteigenvalues     \
         flux_wrappers inviscidfluxlf inviscidfluxzip      \
         inviscidfluxroe inviscidfluxrp viscousfluxcd2     \
         viscousfluxcd4 viscousfluxtcd 

# executables using explicit time-stepping depend on these
EXPSRC = mhdsolverk4 expRMHD_driver

# executables using KINSOL depend on these
KINSRC = impRMHD_interface impRMHD_driver nvector_mhd fnvector_mhd

# executables using CVODE depend on these
CVSRC = impCVODE_interface impCVODE_driver nvector_mhd fnvector_mhd

# executables using PETSc depend on these
PETSCSRC = petsc_interface petsc_driver

# header files (for dependency check when recompiling)
HEADS = nvector_mhd.h mesh.h mesh_common.h mesh_parms.h  \
        mesh_parms_static.h mesh_uparms.h properties.h   \
        mesh_uparms_static.h mgparams.h mpistuff.h boundaries.h
HEADERS = $(addprefix include/, $(HEADS))

# expStatic source files and corresp. object files
SRC1 = $(ALLSRC) $(EXPSRC) init_static
OBJ1 = $(addprefix source/, $(addsuffix .o, $(SRC1)))

# kinStatic source files and corresp object files
SRC2 = $(ALLSRC) $(KINSRC) init_static
OBJ2 = $(addprefix source/, $(addsuffix .o, $(SRC2)))

# cvStatic source files and corresp object files
SRC4 = $(ALLSRC) $(CVSRC) init_static
OBJ4 = $(addprefix source/, $(addsuffix .o, $(SRC4)))

# expLinWave source files and corresp. object files
SRC6 = $(ALLSRC) $(EXPSRC) initwave3
OBJ6 = $(addprefix source/, $(addsuffix .o, $(SRC6)))

# kinLinWave source files and corresp. object files
SRC8 = $(ALLSRC) $(KINSRC) initwave3
OBJ8 = $(addprefix source/, $(addsuffix .o, $(SRC8)))

# cvLinWave source files and corresp. object files
SRC7 = $(ALLSRC) $(CVSRC) initwave3
OBJ7 = $(addprefix source/, $(addsuffix .o, $(SRC7)))

# expKH source files and corresp. object files
SRC9 = $(ALLSRC) $(EXPSRC) initKH
OBJ9 = $(addprefix source/, $(addsuffix .o, $(SRC9)))

# kinKH source files and corresp. object files
SRC11 = $(ALLSRC) $(KINSRC) initKH 
OBJ11 = $(addprefix source/, $(addsuffix .o, $(SRC11)))

# cvKH source files and corresp. object files
SRC10 = $(ALLSRC) $(CVSRC) initKH
OBJ10 = $(addprefix source/, $(addsuffix .o, $(SRC10)))

# exp_island source files and corresp. object files
SRC12 = $(ALLSRC) $(EXPSRC) init_island
OBJ12 = $(addprefix source/, $(addsuffix .o, $(SRC12)))

# kin_island source files and corresp. object files
SRC13 = $(ALLSRC) $(KINSRC) init_island
OBJ13 = $(addprefix source/, $(addsuffix .o, $(SRC13)))

# cv_island source files and corresp. object files
SRC14 = $(ALLSRC) $(CVSRC) init_island
OBJ14 = $(addprefix source/, $(addsuffix .o, $(SRC14)))

# petsc_island source files and corresp. object files
SRC15 = $(ALLSRC) $(PETSCSRC) init_island
OBJ15 = $(addprefix source/, $(addsuffix .o, $(SRC15)))

# fnvec_mhd_test source files and corresp. object files
FNVEC_TEST_SRC = source/modules source/nvector_mhd source/fnvector_mhd \
                 tests/test_nvec_mhd tests/test_fnvec_mhd 
FNVEC_TEST_OBJ = $(addsuffix .o, $(FNVEC_TEST_SRC))



####### Remaining definitions (no adjustment necessary) #######

# SUNDIALS interface requirements
SUND_INCD   = -I${SUNDROOT}/include
SUND_LIBD   = -L${SUNDROOT}/lib
SUND_LIBS   = -lm
CVODE_LIBS  = -lsundials_fcvode -lsundials_cvode 
KINSOL_LIBS = -lsundials_fkinsol -lsundials_kinsol

# Final usable variables
INCS    = -Iinclude ${MPI_INCD} ${SUND_INCD} ${C_INCD} ${PETSC_INCD} ${SILO_INCD}
KLIBS   = ${SUND_LIBD} ${KINSOL_LIBS} ${SUND_LIBS} ${MPI_LIBD} ${MPI_LIBS} ${SILO_LIBS}
CVLIBS  = ${SUND_LIBD} ${CVODE_LIBS} ${SUND_LIBS} ${MPI_LIBD} ${MPI_LIBS} ${SILO_LIBS}
KLIBSH  = ${SUND_LIBD} ${KINSOL_LIBS} ${SUND_LIBS} ${MPI_LIBD} ${MPI_LIBS} ${SILO_LIBS}
CVLIBSH = ${SUND_LIBD} ${CVODE_LIBS} ${SUND_LIBS} ${MPI_LIBD} ${MPI_LIBS} ${SILO_LIBS}

# compilation cleanup files
TRASH = source/*.o tests/*.o *.mod source/*.il



####### Makefile targets #######

help :
	@echo " "
	@echo "usage: 'make <command>', where <command> is one of:"
	@echo " "
	@echo "  Static Tests (resistive or ideal):"
	@echo "     expStatic         [explicit method, fixed dt]"
	@echo "     kinStatic         [KINSOL, fixed dt]"
	@echo "     cvStatic          [CVODE, variable dt]"
	@echo " "
	@echo "  Kelvin-Helmholtz Tests (ideal):"
	@echo "     expKH             [explicit, fixed dt]"
	@echo "     kinKH             [KINSOL, fixed dt]"
	@echo "     cvKH              [CVODE, fixed dt]"
	@echo " "
	@echo "  Linear Wave Tests (ideal):"
	@echo "     expLinWave        [explicit, fixed dt]"
	@echo "     kinLinWave        [KINSOL, fixed dt]"
	@echo "     cvLinWave         [CVODE, variable dt]"
	@echo " "
	@echo " "
	@echo "  Island Coalescence Tests (resistive):"
	@echo "     exp_island        [explicit, fixed dt]"
	@echo "     kin_island        [KINSOL, fixed dt]"
	@echo "     cv_island         [CVODE, variable dt]"
	@echo "     petsc_island      [PETSc]"
	@echo " "
	@echo "  Code Tests:"
	@echo "     fnvec_mhd_test    [checks F90/C vector kernel interface]"
	@echo " "
	@echo "  Clean-up:"
	@echo "     clean             [cleans temporary build files]"
	@echo "     outclean          [cleans output from a run]"
	@echo "     tempclean         [cleans temporary build and ~ files]"
	@echo " "
	@echo "Note: controls over spatial dimension, parallel execution,"
	@echo "      boundary conditions, dynamic meshing, upwind vs CD"
	@echo "      inviscid fluxes, order of CD inviscid fluxes,"
	@echo "      primitive vs conservative upwind fluxing, use of"
	@echo "      Harten's entropy fix, and use of viscous fluxes"
	@echo "      are all decided through compiler definition flags"
	@echo "      in the file Make.config.  See that file for info."

expStatic : ${HEADERS} ${OBJ1}
	${F90} ${FFLAGS_} -o $@ ${OBJ1} ${INCS} ${MPI_LIBD} ${MPI_LIBS} 

kinStatic : ${HEADERS} ${OBJ2}
	${F90} ${FFLAGS_} -o $@ ${OBJ2} ${INCS} ${KLIBS} 

kinStatic2 : ${HEADERS} ${OBJ2N}
	${F90} ${FFLAGS_} -o $@ ${OBJ2N} ${INCS} ${KLIBS} 

cvStatic : ${HEADERS} ${OBJ4}
	${F90} ${FFLAGS_} -o $@ ${OBJ4} ${INCS} ${CVLIBS} 

cvStatic2 : ${HEADERS} ${OBJ4N}
	${F90} ${FFLAGS_} -o $@ ${OBJ4N} ${INCS} ${CVLIBS} 

expLinWave : ${HEADERS} ${OBJ6}
	${F90} ${FFLAGS_} -o $@ ${OBJ6} ${INCS} ${MPI_LIBD} ${MPI_LIBS} 

kinLinWave : ${HEADERS} ${OBJ8}
	${F90} ${FFLAGS_} -o $@ ${OBJ8} ${INCS} ${KLIBS} 

kinLinWave2 : ${HEADERS} ${OBJ8N}
	${F90} ${FFLAGS_} -o $@ ${OBJ8N} ${INCS} ${KLIBS} 

cvLinWave : ${HEADERS} ${OBJ7}
	${F90} ${FFLAGS_} -o $@ ${OBJ7} ${INCS} ${CVLIBS} 

cvLinWave2 : ${HEADERS} ${OBJ7N}
	${F90} ${FFLAGS_} -o $@ ${OBJ7N} ${INCS} ${CVLIBS} 

expKH : ${HEADERS} ${OBJ9}
	${F90} ${FFLAGS_} -o $@ ${OBJ9} ${INCS} ${MPI_LIBD} ${MPI_LIBS} 

kinKH : ${HEADERS} ${OBJ11}
	${F90} ${FFLAGS_} -o $@ ${OBJ11} ${INCS} ${KLIBS} 

cvKH : ${HEADERS} ${OBJ10}
	${F90} ${FFLAGS_} -o $@ ${OBJ10} ${INCS} ${CVLIBS} 

exp_island : ${HEADERS} ${OBJ12}
	${F90} ${FFLAGS_} -o $@ ${OBJ12} ${INCS} ${MPI_LIBD} ${MPI_LIBS} 

kin_island : ${HEADERS} ${OBJ13}
	${F90} ${FFLAGS_} -o $@ ${OBJ13} ${INCS} ${KLIBS} 

cv_island : ${HEADERS} ${OBJ14}
	${F90} ${FFLAGS_} -o $@ ${OBJ14} ${INCS} ${CVLIBS} 

petsc_island : ${HEADERS} ${OBJ15}
	${F90} ${FFLAGS_} -o $@ ${OBJ15} ${INCS} ${PETSC_LIBS}

fnvec_mhd_test : ${HEADERS} ${FNVEC_TEST_OBJ}
	${F90} ${FFLAGS_} -o $@ ${FNVEC_TEST_OBJ} ${INCS} ${KLIBS}


clean:
	\rm -f ${TRASH}

outclean:
	\rm -rf *.txt *.history data* dump.* gmon.out *.png *.grid.* multimesh* dump*.bin

tempclean: clean
	\rm -i *~ source/*~ include/*~



####### Makefile pattern rule redefinitions #######

.F.o :
	${F90} -c ${FFLAGS_} ${MPI_INCD} ${INCS} $< -o $@

.f.o :
	${F77} -c ${F77FLAGS} ${INCS} $< -o $@

.c.o :
	${CC} -c ${CFLAGS_} ${INCS} $< -o $@


####### End of Makefile #######
