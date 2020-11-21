CC   = gcc -Wall -O3
OBJS = main.o pareto_set.o bounds.o
PROG = prog

# GENERIC

all:	 ${PROG}

clean:
		 rm ${OBJS} ${PROG}

${PROG}: ${OBJS}
		 ${CC} ${OBJS} -o $@
.c.o:
		 ${CC} $< -c -o $@

#############################

pareto_set.o:	pareto_set.c pareto_set.h
bounds.o:	bounds.c bounds.h
main.o: main.c pareto_set.h bounds.h
