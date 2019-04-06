#local
export FC = gfortran
#export FFLAGS = -O3 
export FFLAGS = -g -fcheck=bounds -Wall -Wextra -Wshadow
export bindir = ../bin

objects = input.o vcount.o levels.o recur.o sort.o

all: $(objects)
	if [ ! -d $(bindir) ]; then mkdir -p $(bindir); fi
	cp vcount.o $(bindir)/vcount

input.o : input.f90 
	$(FC) $(FFLAGS) -c input.f90

levels.o : levels.f90 recur.o
	$(FC) $(FFLAGS) -c levels.f90 recur.o

recur.o : recur.f90
	$(FC) $(FFLAGS) -c recur.f90

sort.o : sort.f90
	$(FC) $(FFLAGS) -c sort.f90

vcount.o : vcount.f90 input.o levels.o recur.o sort.o
	$(FC) $(FFLAGS) -o vcount.o vcount.f90 input.o levels.o recur.o sort.o

clean:
	rm *.o *.mod


