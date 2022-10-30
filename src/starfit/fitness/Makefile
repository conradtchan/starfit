SOURCE ?= .

LIBRARY_OBJECTS = \
	mleqs.o \
	norm.o \
	powell.o \
	solver.o \
	type_def.o

F90 = gfortran -c -Ofast -fno-finite-math-only -fPIC -Wall

a.out	: ${LIBRARY_OBJECTS}
	gfortran ${LIBRARY_OBJECTS} -o a.out

type_def.o type_def.mod &: ${SOURCE}/type_def.f90
	@rm -f type_def.mod
	$(F90) $(FFLAGS) $<

norm.o norm.mod &: ${SOURCE}/norm.f90 \
	type_def.mod
	@rm -f morm.mod
	$(F90) $(FFLAGS) $<

powell.o powell.mod &: ${SOURCE}/powell.f90 \
	type_def.mod
	@rm -f morm.mod
	$(F90) $(FFLAGS) $<

mleqs.o mleqs.mod &: ${SOURCE}/mleqs.f90 \
	type_def.mod
	@rm -f mleqs.mod
	$(F90) $(FFLAGS) $<

solver.o solver.mod star_data.mod abu_data.mod &: ${SOURCE}/solver.f90 \
	type_def.mod
	@rm -f solver.mod abu_data.mod star_data.mod
	$(F90) $(FFLAGS) $<

.PHONY:	clean

clean:
	-rm -f *.o *.a *.mod *.smod *~ \#*\# .*~ .\#*
	-rm -f a.out
