#CXX=gcc
#FC=gfortran
#CXX=icpc
#FC=ifort -std08

include $(KOKKOS_ROOT)/Makefile.kokkos

KOKKOS_INC=$(KOKKOS_CXXFLAGS)
KOKKOS_LIB=$(KOKKOS_CXXLDFLAGS) $(KOKKOS_LIBS) $(KOKKOS_LINK_FLAGS)

OMP = -fopenmp
DEBUG = -g
SRCDIR=../src
LIBFLCLROOT=../lib/flcl
LIBFLCLBUILD=$(LIBFLCLROOT)/build/

ifeq ($(FC),ifort)
  FSTD = -std08
  OMP = -qopenmp
endif
ifeq ($(FC),gfortran)
  FSTD = -std=f2008
  OMP = -fopenmp
endif
ifeq ($(KOKKOS_DEVICES),Cuda)
  CXX=$(NVCC_WRAPPER)
endif

all: xkt.x

amr-f.o: $(SRCDIR)/interface/amr-f.f90
	$(FC) $(FSTD) $(DEBUG) -I$(LIBFLCLBUILD) -I$(SRCDIR)/interface -c $(SRCDIR)/interface/amr-f.f90

interface-f.o: $(SRCDIR)/interface/interface-f.f90
	$(FC) $(FSTD) $(DEBUG) -I$(LIBFLCLBUILD) -I$(SRCDIR)/interface -c $(SRCDIR)/interface/interface-f.f90

interface-cxx.o: $(SRCDIR)/interface/interface-cxx.cpp
	$(CXX) $(DEBUG) -I$(LIBFLCLROOT)/src -I$(SRCDIR)/levels $(KOKKOS_INC) -c $(SRCDIR)/interface/interface-cxx.cpp

levels-f.o: $(SRCDIR)/levels/levels-f.f90
	$(FC) $(FSTD) $(DEBUG) -I$(LIBFLCLBUILD) -I$(SRCDIR)/levels -c $(SRCDIR)/levels/levels-f.f90

main.o: $(SRCDIR)/main/main.f90
	$(FC) $(FSTD) $(DEBUG) -I$(LIBFLCLBUILD) -c $(SRCDIR)/main/main.f90

xkt.x: levels-f.o  interface-f.o amr-f.o interface-cxx.o main.o
	# $(FC)  levels.o abi.o abi-cxx.o f_interface.o f_amr.o c_interface.o main.o $(KOKKOS_LIB) -lstdc++ $(OMP) -o xkt.x
	$(FC) $(FSTD) levels-f.o interface-f.o amr-f.o interface-cxx.o main.o $(KOKKOS_LIB) -L$(LIBFLCLBUILD) -lflcl -lstdc++ -o xkt.x


clean:
	rm -f ./*.o ./*.mod ./xkt.x
