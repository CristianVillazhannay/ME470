ALL: ex2
CLEANFILES: ex2

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

ex2: ex2.o chkopts
	${CLINKER} -o ex2 ex2.o ${PETSC_LIB}
	${RM} ex2.o

