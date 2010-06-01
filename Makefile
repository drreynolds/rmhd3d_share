###################################################################
#  Makefile for unigrid resistive MHD code
#  Daniel R. Reynolds, UCSD Mathematics
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
FFLAGS_ = $(FFLAGS) $(PROB_DEFS) $(SUND_DEFS)
CFLAGS_ = $(CFLAGS) $(PROB_DEFS) $(SUND_DEFS)



####### Executable source dependencies #######

# all executables depend on these files
ALLSRC = modules profiling main setupgrid setuplocaldomain \
         setboundaryvalues output bdryexchange             \
         seteigenvectors_cons inviscidfluxcd2              \
         inviscidfluxcd4 inviscidfluxtcd projection newdt  \
         xchange diags drivertools mhdtools derf

# executables including viscous fluxes depend on these
VISCSRC = viscousfluxcd2 viscousfluxcd4 viscousfluxtcd 

# executables using explicit time-stepping depend on these
EXPSRC = mhdsolverk4 expRMHD_driver

# executables using KINSOL depend on these
KINSRC = impRMHD_interface impRMHD_driver nvector_mhd fnvector_mhd

# executables using CVODE depend on these
CVSRC = impCVODE_interface impCVODE_driver nvector_mhd fnvector_mhd

# header files (for dependency check when recompiling)
HEADS = nvector_mhd.h mesh.h mesh_common.h mesh_parms.h  \
        mesh_parms_static.h mesh_uparms.h properties.h   \
        mesh_uparms_static.h mgparams.h mpistuff.h boundaries.h
HEADERS = $(addprefix include/, $(HEADS))

# necessary lapack/blas files
SRC_LPK = dgemm dgemv dger dgetf2 dgetrf dgetrs dlaswp dscal \
          dswap dtrsm idamax ieeeck ilaenv xerbla lsame
OBJ_LPK = $(addprefix source/Lapack/, $(addsuffix .o, $(SRC_LPK)))

# preconditioning files (hyperbolic and diffusive, respectively)
PRECSRCH = ptri_parallel fastwave_prec combo_prec
PRECSRCD = viscous_prec-Du viscous_prec-Db viscous_prec-De \
           vprec_solver vprec_mult combo_prec

# expStatic source files and corresp. object files
SRC1 = $(ALLSRC) $(VISCSRC) $(EXPSRC) init_static
OBJ1 = $(addprefix source/, $(addsuffix .o, $(SRC1)))

# kinStatic source files and corresp object files
SRC2 = $(ALLSRC) $(VISCSRC) $(PRECSRCH) $(PRECSRCD) $(KINSRC) init_static
OBJ2 = $(addprefix source/, $(addsuffix .o, $(SRC2))) $(OBJ_LPK)

# cvStatic source files and corresp object files
SRC4 = $(ALLSRC) $(VISCSRC) $(PRECSRCH) $(PRECSRCD) $(CVSRC) init_static
OBJ4 = $(addprefix source/, $(addsuffix .o, $(SRC4))) $(OBJ_LPK)

# expLinWave source files and corresp. object files
SRC6 = $(ALLSRC) $(EXPSRC) initwave3
OBJ6 = $(addprefix source/, $(addsuffix .o, $(SRC6)))

# kinLinWave source files and corresp. object files
SRC8 = $(ALLSRC) $(PRECSRCH) $(KINSRC) initwave3
OBJ8 = $(addprefix source/, $(addsuffix .o, $(SRC8))) $(OBJ_LPK)

# cvLinWave source files and corresp. object files
SRC7 = $(ALLSRC) $(PRECSRCH) $(CVSRC) initwave3
OBJ7 = $(addprefix source/, $(addsuffix .o, $(SRC7))) $(OBJ_LPK)

# expKH source files and corresp. object files
SRC9 = $(ALLSRC) $(VISCSRC) $(EXPSRC) initKH
OBJ9 = $(addprefix source/, $(addsuffix .o, $(SRC9)))

# kinKH source files and corresp. object files
SRC11 = $(ALLSRC) $(VISCSRC) $(PRECSRCH) $(PRECSRCD) $(KINSRC) initKH 
OBJ11 = $(addprefix source/, $(addsuffix .o, $(SRC11))) $(OBJ_LPK)

# cvKH source files and corresp. object files
SRC10 = $(ALLSRC) $(VISCSRC) $(PRECSRCH) $(PRECSRCD) $(CVSRC) initKH
OBJ10 = $(addprefix source/, $(addsuffix .o, $(SRC10))) $(OBJ_LPK)

# test_driver source files and corresp. object files
TEST_DRIVER_SRC = $(ALLSRC) $(PRECSRCH) $(PRECSRCD) init2 test_driver \
                  impCVODE_interface nvector_mhd fnvector_mhd
TEST_DRIVER_OBJ = $(addprefix source/, $(addsuffix .o, $(TEST_DRIVER_SRC))) $(OBJ_LPK)

# fnvec_mhd_test source files and corresp. object files
FNVEC_TEST_SRC = source/modules source/nvector_mhd source/fnvector_mhd \
                 tests/test_nvec_mhd tests/test_fnvec_mhd 
FNVEC_TEST_OBJ = $(addsuffix .o, $(FNVEC_TEST_SRC))

# hypre_test source files and corresp. object files
HYPRE_TEST_SRC = hypre_test viscous_prec-Du    \
                 viscous_prec-Db viscous_prec-De
HYPRE_TEST_OBJ = $(addprefix source/, $(addsuffix .o, $(HYPRE_TEST_SRC)))

# fhypre_test source files and corresp. object files
FHYPRE_TEST_SRC = fhypre_test viscous_prec-Du    \
                  viscous_prec-Db viscous_prec-De
FHYPRE_TEST_OBJ = $(addprefix source/, $(addsuffix .o, $(FHYPRE_TEST_SRC)))



####### Remaining definitions (no adjustment necessary) #######

# SUNDIALS interface requirements
SUND_INCD   = -I${SUNDROOT}/include
SUND_LIBD   = -L${SUNDROOT}/lib
SUND_LIBS   = -lm
CVODE_LIBS  = -lsundials_fcvode -lsundials_cvode 
KINSOL_LIBS = -lsundials_fkinsol -lsundials_kinsol

# HYPRE interface requirements
HYPRE_INCD   = -I${HYPREROOT}/include
HYPRE_LIBD   = -L${HYPREROOT}/lib
HYPRE_LIBS   = -lHYPRE

# Final usable variables
INCS    = -Iinclude ${MPI_INCD} ${SUND_INCD} ${HYPRE_INCD} ${C_INCD}
KLIBS   = ${SUND_LIBD} ${KINSOL_LIBS} ${SUND_LIBS} ${HYPRE_LIBD} ${HYPRE_LIBS} ${MPI_LIBD} ${MPI_LIBS}
CVLIBS  = ${SUND_LIBD} ${CVODE_LIBS} ${SUND_LIBS} ${HYPRE_LIBD} ${HYPRE_LIBS} ${MPI_LIBD} ${MPI_LIBS}
KLIBSH   = ${SUND_LIBD} ${KINSOL_LIBS} ${SUND_LIBS} ${MPI_LIBD} ${MPI_LIBS}
CVLIBSH  = ${SUND_LIBD} ${CVODE_LIBS} ${SUND_LIBS} ${MPI_LIBD} ${MPI_LIBS}

# compilation cleanup files
TRASH = source/*.o source/Lapack/*.o tests/*.o *.mod source/*.il



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
	@echo "  Code Tests:"
	@echo "     test_driver       [driver to run customized tests]"
	@echo "     fnvec_mhd_test    [checks F90/C vector kernel interface]"
	@echo "     hypre_test        [checks C hypre interface]"
	@echo "     fhypre_test       [checks F90 hypre interface]"
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

kinKH2 : ${HEADERS} ${OBJ11N}
	${F90} ${FFLAGS_} -o $@ ${OBJ11N} ${INCS} ${KLIBS} 

cvKH : ${HEADERS} ${OBJ10}
	${F90} ${FFLAGS_} -o $@ ${OBJ10} ${INCS} ${CVLIBS} 

cvKH2 : ${HEADERS} ${OBJ10N}
	${F90} ${FFLAGS_} -o $@ ${OBJ10N} ${INCS} ${CVLIBS} 

test_driver : ${HEADERS} ${TEST_DRIVER_OBJ}
	${F90} ${FFLAGS_} -o $@ ${TEST_DRIVER_OBJ} ${INCS} ${CVLIBS} 

fnvec_mhd_test : ${HEADERS} ${FNVEC_TEST_OBJ}
	${F90} ${FFLAGS_} -o $@ ${FNVEC_TEST_OBJ} ${INCS} ${KLIBS}

hypre_test : ${HEADERS} ${HYPRE_TEST_OBJ}
	${F90} ${FFLAGS_} -o $@ ${HYPRE_TEST_OBJ} ${INCS} ${HYPRE_LIBD} ${HYPRE_LIBS} ${MPI_LIBD} ${MPI_LIBS} -lm

fhypre_test : ${HEADERS} ${FHYPRE_TEST_OBJ}
	${F90} ${FFLAGS_} -o $@ ${FHYPRE_TEST_OBJ} ${INCS} ${HYPRE_LIBD} ${HYPRE_LIBS} ${MPI_LIBD} ${MPI_LIBS} -lm


clean:
	\rm -f ${TRASH}

outclean:
	\rm -f C_test* F_test* *.txt *.history output.0* outavs.00* dump.* gmon.out precmat.* rhsvec.*

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
