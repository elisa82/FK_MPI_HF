##  If you want to change the size of GF bank, please just modify the size.h
F77=mpif90 -O3 
#F77=mpif90 -O3 -fallow-argument-mismatch
#F77=/usr/lib64/mpich/bin/mpif90 -O3 -f90=ifort -assume byterecl
CC =gcc -O3

SUBS = fft.o Complex.o distaz.o 
FKSUBS = gfbank_mpi.o cmodel.o fk.o kernel.o compound.o haskell.o source.o bessel.o $(SUBS)
all: gfbank_mpi clean

gfbank_mpi: \
	$(FKSUBS)
	$(F77) -o gfbank_mpi $(FKSUBS)
fft.o:\
	fft.c 
	$(CC)  -c fft.c

Complex.o:\
	Complex.c Complex.h
	$(CC)  -c Complex.c

distaz.o:\
	distaz.f
	$(F77) -c  distaz.f

cmodel.o:\
	cmodel.f
	$(F77) -c  cmodel.f

fk.o:\
	fk.f
	$(F77) -c fk.f

kernel.o:\
	kernel.f
	$(F77) -c kernel.f

compound.o:\
	compound.f
	$(F77) -c compound.f

haskell.o:\
	haskell.f
	$(F77) -c haskell.f

source.o:\
	source.f
	$(F77) -c source.f
bessel.o:\
	bessel.f
	$(F77) -c bessel.f

gfbank_mpi.o: \
	gfbank_mpi.f
	$(F77) -c gfbank_mpi.f
syn_fk: syn_fk.f sacio.o
	$(F77) -o $@ $@.f sacio.o

clean:
	rm -f *.o
