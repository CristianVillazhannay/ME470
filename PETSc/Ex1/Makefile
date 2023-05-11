ALL: ex1
CLEANFILES = ex1

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

ex1: ex1.o chkopts
	${CLINKER} -o ex1 ex1.o ${PETSC_LIB}
	${RM} ex1.o

