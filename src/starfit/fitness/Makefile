SOURCE ?= .

LIBRARY_OBJECTS = \
	type_def.o \
	mleqs.o \
	norm.o \
	powell.o \
	abu_data.o \
	star_data.o \
	solver.o \


FFLAGS = -c -Ofast -fno-finite-math-only -fPIC -Wall

F90 = gfortran

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
	@rm -f powell.mod
	$(F90) $(FFLAGS) $<

mleqs.o mleqs.mod &: ${SOURCE}/mleqs.f90 \
	type_def.mod
	@rm -f mleqs.mod
	$(F90) $(FFLAGS) $<

abu_data.o abu_data.mod &: ${SOURCE}/abu_data.f90 \
	type_def.mod
	@rm -f abu_data.mod
	$(F90) $(FFLAGS) $<

star_data.o star_data.mod &: ${SOURCE}/star_data.f90 \
	type_def.mod mleqs.mod
	@rm -f star_data.mod
	$(F90) $(FFLAGS) $<

solver.o solver.mod &: ${SOURCE}/solver.f90 \
	type_def.mod powell.mod mleqs.mod norm.mod \
	star_data.mod abu_data.mod
	@rm -f solver.mod
	$(F90) $(FFLAGS) $<

.PHONY:	clean

clean:
	-rm -f *.o *.a *.mod *.smod *~ \#*\# .*~ .\#*
	-rm -f a.out
