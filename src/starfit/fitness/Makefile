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
	typedef.o \
	utils.o \
	linalg.o \
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

typedef.o typedef.mod &: ${SOURCE}/typedef.f90
	@rm -f typedef.mod
	$(F90) $(FFLAGS) $<

utils.o utils.mod &: ${SOURCE}/utils.f90 \
	typedef.mod
	@rm -f utils.mod
	$(F90) $(FFLAGS) $<

norm.o norm.mod &: ${SOURCE}/norm.f90 \
	typedef.mod
	@rm -f morm.mod
	$(F90) $(FFLAGS) $<

powell.o powell.mod &: ${SOURCE}/powell.f90 \
	typedef.mod
	@rm -f powell.mod
	$(F90) $(FFLAGS) $<

linalg.o linalg.mod &: ${SOURCE}/linalg.f90 \
	typedef.mod utils.mod
	@rm -f linalg.mod
	$(F90) $(FFLAGS) $<

abu_data.o abu_data.mod &: ${SOURCE}/abu_data.f90 \
	typedef.mod
	@rm -f abu_data.mod
	$(F90) $(FFLAGS) $<

star_data.o star_data.mod &: ${SOURCE}/star_data.f90 \
	typedef.mod linalg.mod utils.mod
	@rm -f star_data.mod
	$(F90) $(FFLAGS) $<

fitting.o fitting.mod &: ${SOURCE}/fitting.f90 \
	typedef.mod powell.mod linalg.mod norm.mod \
	star_data.mod abu_data.mod utils.mod
	@rm -f fitting.mod
	$(F90) $(FFLAGS) $<

tests.o tests.mod &: ${SOURCE}/tests.f90 \
	typedef.mod fitting.mod
	@rm -f tests.mod
	$(F90) $(FFLAGS) $<

test.o &: ${SOURCE}/test.f90 \
	typedef.mod tests.mod
	$(F90) $(FFLAGS) $<

test_leqs.o &: ${SOURCE}/test_leqs.f90 \
	typedef.mod linalg.mod
	$(F90) $(FFLAGS) $<

test_leqs: test_leqs.o linalg.o utils.o
	${LINK} ${LFLASG} test_leqs.o linalg.o utils.o -o $@

.PHONY:	clean, ar

clean:
	-rm -f *.o *.a *.mod *.smod *~ \#*\# .*~ .\#* _*.f _*.c fortranobject.c fortranobject.h _*.f90
	-rm -f ${EXECUTABLE} ${LIBRARY}

${LIBRARY}: ${LIBRARY_OBJECTS}
	${AR} ${AFLAGS} ${LIBRARY} ${LIBRARY_OBJECTS}

ar: ${LIBRARY}
