%.o: %.mod
MF=Makefile

# compiler and flags
#------------------------------------------------------------------------

EXE=elast.x
# GNU compiler
FC=gfortran
CFLAGS = -fallow-argument-mismatch -O3 $(UFLAG) $(DEBUG) 
# intel compiler
#FC=ifort -mkl=sequential

# use your own blas and lapack libs
BLAS=-L../lib/BLAS-3.10.0 -lblas
LAPACK=-L../lib/lapack-3.10.0 -llapack 
LPKINC=-I../lib/lapack-3.10.0/LAPACKE/include/

# or use system compiler environment (GNUi/intel fortran compiler)
#BLAS=
#LAPACK=
#LPKINC=

# or on Mac OS, use Apple's Accelerate
#LAPACK=  -framework Accelerate

# User flag
UFLAG =
# Debug flags
DEBUG = 
#-----------------------------------------------------------------------
# do not modify the following
#-----------------------------------------------------------------------

SRC= \
	nrtype.f90 data.f90 dlexlib.f90 doread.f90 readin.f90 dist.f90 \
	dev.f90 extract.f90 calc.f90 ssnpt.f90 fitting.f90 elast.f90 

.SUFFIXES: .o .f90

OBJ= $(SRC:.f90=.o)

all: $(EXE) $(OBJ)

$(EXE): $(OBJ)
	$(FC) $(OBJ) $(LAPACK) $(LPKLINC) $(BLAS) -o $@ $(CFLAGS)
	@echo "\o/ Compilation was successful!"

.f90.o:
	$(FC) $(LAPACK) $(BLAS) $(BLAS) -c $< $(CFLAGS)

$(OBJ): $(MF) 

clean:
	rm *.o *.mod $(EXE)
