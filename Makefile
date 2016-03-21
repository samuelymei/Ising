.SUFFIXES: 
.SUFFIXES: .f90 .o

FC = mpif90
FFLAGS = -CB
INCLUDE = 
LIBS = 

EXE = Ising.x

MODULES = constant.mod precision_m.mod random_m.mod ising2d_m.mod

OBJS = precision_m.o lib.o constant.o random.o ising2d.o Ising.o

all:	${EXE}

$(EXE):$(OBJS) ${MODULES}
	$(FC) -o $@ $(FFLAGS) $(OBJS) $(LIBS)

%.o %.mod:%.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $<

include .depend

depend .depend:
	makedepf90 *.f90 > .depend

clean:
	/bin/rm -f $(EXE) $(OBJS) ${MODULES}

