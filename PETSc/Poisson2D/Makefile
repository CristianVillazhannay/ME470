ALL: petscPoisson2D
CLEANFILES = petscPoisson2D

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

petscPoisson2D: petscPoisson2D.o chkopts
	${CLINKER} -o petscPoisson2D petscPoisson2D.o ${PETSC_LIB}
	${RM} petscPoisson2D.o
