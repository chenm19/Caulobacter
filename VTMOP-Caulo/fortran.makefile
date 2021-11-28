# define matlab dir
MDIR = /Applications/MATLAB_R2020a.app
# for linux
#MDIR = /etc/matlab 

# compiles mex files using gfortran
F90 = gfortran
F77 = gfortran
COMP = -c

# Set the compiler and flags here
# F90 = gfortran -std=f2008 -fopenmp
# F77 = gfortran -std=legacy -fopenmp
# Set the build flags here
# COMP = -c
# LIBS = -llapack -lblas

# Add compiler specific flags.
ifeq ($(F90),g90)
	FFLAGS+=-O3 -Wall -ftrace=full -std=f90
else ifeq ($(F90),gfortran)
	FFLAGS+= -O3 -cpp 
# 	FFLAGS+=-O3 -Wall 
	FFLAGS+=-ffpe-trap=invalid,zero,overflow
# 	FFLAGS+=-ffpe-summary=all
# 	FFLAGS+=-fcheck=all -fbacktrace
# Use the following line to enable OpenMP.
# 	FFLAGS+=-fopenmp
# For debugging, uncomment the following line.
# 	FFLAGS+=-g
else ifeq ($(F90),lf90)
	FFLAGS+=-Wa,--32 --f90 --in --trap
else ifeq ($(F90),f90)
	FFLAGS+=-en -O2 -s
	FFLAGS+=-TENV:simd_imask=off -TENV:simd_zmask=off -TENV:simd_omask=off
endif

# to use the intel compiler instead, uncomment F95 and FFLAGS below:

# compiles mex file using the intel compiler
# F95 = ifort

# compiler flags for intel compiler
# FFLAGS = -O3 -fpp -fPIC -D__amd64 -module fortran

# Figure out which platform we're on
UNAME = $(shell uname -s)

# Linux
ifeq ($(findstring Linux,${UNAME}), Linux)
	# define which files to be included
	FINCLUDE = -I$(MDIR)/extern/include -Ifortran -shared 
	# define extension
	EXT = mexa64
endif

# Mac OS X
ifeq ($(findstring Darwin,${UNAME}), Darwin)
	# define which files to be included
	FINCLUDE = -L$(MDIR)/bin/maci64 -Ifortran -lmx -lmex -lmat 
# 	FINCLUDE = -L$(MDIR)/bin/maci64 -I$(MDIR)/extern/include -Ifortran -lmx -lmex -lmat
	# define extension
	EXT = mexmaci64
    ifeq ($(F90), ifort)
        FNOMAIN = -nofor_main -bundle
    endif
    ifeq ($(F90), gfortran)
        FNOMAIN = -shared -fopenmp
    endif
endif



#LIBS = -llapack -lblas
# List of object files
OBJS = vtmop.o vtmop_func.o linear_shepard.o svtdirect.o bvtdirect.o \
	shared_modules.o delsparse.o qnstop.o slatec.o lapack.o blas.o \
	model.o main.o gatewayFortran.o

# Execute the entire project
all : gatewayFortran.$(EXT)  
# 	./gatewayFortran
# 	export OMP_NESTED=TRUE && export OMP_NUM_THREADS=2,2 && ./gatewayFortran

gatewayFortran.$(EXT) : $(OBJS)
	$(F90) $(FFLAGS) $(FINCLUDE) $(FNOMAIN) $(LIBS) -o $@ $^

gatewayFortran.o : gatewayFortran.f90 main.o
	$(F90) $(FFLAGS) $(COMP) gatewayFortran.f90 -o gatewayFortran.o

main.o : main.f90 model.o 
	$(F90) $(FFLAGS) $(COMP) main.f90 -o main.o

model.o : model.f90 
	$(F90) $(FFLAGS) $(COMP) model.f90 -o model.o

vtmop.o : vtmop.f90 delsparse.o linear_shepard.o svtdirect.o bvtdirect.o \
	  qnstop.o
	$(F90) $(FFLAGS) $(COMP) vtmop.f90 -o vtmop.o

vtmop_func.o: vtmop_func.f90 shared_modules.o
	$(F90) $(FFLAGS) $(COMP) vtmop_func.f90 -o vtmop_func.o

delsparse.o : delsparse.f90 shared_modules.o
	$(F90) $(FFLAGS) $(COMP) delsparse.f90 -o delsparse.o

linear_shepard.o : linear_shepard.f90 shared_modules.o
	$(F90) $(FFLAGS) $(COMP) linear_shepard.f90 -o linear_shepard.o

qnstop.o : qnstop.f90 shared_modules.o
	$(F90) $(FFLAGS) $(COMP) qnstop.f90 -o qnstop.o

svtdirect.o : sVTdirect.f90 shared_modules.o
	$(F90) $(FFLAGS) $(COMP) sVTdirect.f90 -o svtdirect.o

bvtdirect.o : bVTdirect.f90 shared_modules.o
	$(F90) $(COMP) bVTdirect.f90 -o bvtdirect.o

shared_modules.o : shared_modules.f90
	$(F90) $(COMP) shared_modules.f90 -o shared_modules.o

slatec.o : slatec.f
	$(F77) $(COMP) slatec.f -o slatec.o

lapack.o : lapack.f
	$(F77) $(COMP) lapack.f -o lapack.o

blas.o : blas.f
	$(F77) $(COMP) blas.f -o blas.o


purge:	clean

clean:
	rm -f $(EXE) *.o *.mod
	rm -f matlab.out 
	rm -f vtmop.dat vtmop.chkpt
	rm -f fort.6 

