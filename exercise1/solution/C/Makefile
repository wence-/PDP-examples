# Leave alone below here


#
# lists of sources, etc
#
SRC = mandelbrot.c

OBJ=$(SRC:.c=.o)

CC = mpicc
CFLAGS=	-Wall -O2 -Iutils

LIBS = utils/arralloc.o utils/write_ppm.o utils/read_options.o

all: mandelbrot

$(OBJ):	$(MAKEFILE) $(SRC) $(INC)

utils:
	$(MAKE) -C utils

mandelbrot:	$(OBJ) utils
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LIBS)


clean: clean-utils
	rm -f mandelbrot *~ $(OBJ) $(EXE) core

clean-utils:
	$(MAKE) -C utils clean

.PHONY: utils clean-utils
