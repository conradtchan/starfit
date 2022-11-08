SOURCE ?= .

#FFLAGS ?= -c -Ofast -fno-finite-math-only -fPIC -Wall
#FFLAGS ?= -c -Ofast -fno-finite-math-only -fno-second-underscore -fPIC -Wall

FFLAGS ?= -c -g -Wall

F90 ?= gfortran

LINK ?= gfortran

LFLAGS ?=

AR ?= ar

AFLAGS ?= crv

LIBRARY = solver.a

LIBRARY_OBJECTS = \
	type_def.o \
	utils.o \
	mleqs.o \
	norm.o \
	powell.o \
	abu_data.o \
	star_data.o \
	fitting.o

EXECUTABLE_OBJECTS = \
	test.o \
	tests.o

EXECUTABLE = \
	test

#${EXECUTABLE}: ${EXECUTABLE_OBJECTS} ${LIBRARY}
#	${LINK} ${LFLASG} -o ${EXECUTABLE} ${EXECUTABLE_OBJECTS} ${LIBRARY}

${EXECUTABLE}: ${EXECUTABLE_OBJECTS} ${LIBRARY_OBJECTS}
	${LINK} ${LFLASG} ${LIBRARY_OBJECTS} ${EXECUTABLE_OBJECTS} -o ${EXECUTABLE}

type_def.o type_def.mod &: ${SOURCE}/type_def.f90
	@rm -f type_def.mod
	$(F90) $(FFLAGS) $<

utils.o utils.mod &: ${SOURCE}/utils.f90 \
	type_def.mod
	@rm -f utils.mod
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
	type_def.mod mleqs.mod utils.mod
	@rm -f star_data.mod
	$(F90) $(FFLAGS) $<

fitting.o fitting.mod &: ${SOURCE}/fitting.f90 \
	type_def.mod powell.mod mleqs.mod norm.mod \
	star_data.mod abu_data.mod utils.mod
	@rm -f fitting.mod
	$(F90) $(FFLAGS) $<

# solver.o solver.mod &: ${SOURCE}/solver.f90 \
# 	type_def.mod fitting.mod
# 	@rm -f solver.mod
# 	$(F90) $(FFLAGS) $<

tests.o tests.mod &: ${SOURCE}/tests.f90 \
	type_def.mod fitting.mod
	@rm -f tests.mod
	$(F90) $(FFLAGS) $<

test.o &: ${SOURCE}/test.f90 \
	type_def.mod tests.mod
	$(F90) $(FFLAGS) $<

.PHONY:	clean, ar

clean:
	-rm -f *.o *.a *.mod *.smod *~ \#*\# .*~ .\#* _*.f _*.c fortranobject.c fortranobject.h _*.f90
	-rm -f ${EXECUTABLE} ${LIBRARY}

${LIBRARY}: ${LIBRARY_OBJECTS}
	${AR} ${AFLAGS} ${LIBRARY} ${LIBRARY_OBJECTS}

ar: ${LIBRARY}
