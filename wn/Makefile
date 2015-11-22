FC=ifort
MPIF90=mpif90
FFLAGS=-I. -O3 -mcmodel=large -shared-intel -w -mkl -fpp -extend-source -finline-functions -no-inline-factor -qopenmp #-ipo
MPIFFLAGS=$(FFLAGS) -DUSE_MPI
#LDFLAGS=-L. -llapack -lblas
SOURCES=$(wildcard *.f90)
OBJECTS=$(SOURCES:.f90=.o)
TARGET=wn

.PHONY: all

all: $(TARGET)

mpi: FC=$(MPIF90) all

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

clean:
	rm *.o *.mod

%.o: %.f90
	$(FC) $(FFLAGS) -c $(@:.o=.f90)

# Dependencies of files
main.o: \
    main.f90

