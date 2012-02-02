# Leave alone below here


#
# lists of sources, etc
#
SRC = fractal-set.c

EXE = fractal-set
OBJ=$(SRC:.c=.o)

include ../../../config.mak

LIBDIR = ../../utils
CFLAGS := $(CFLAGS) -I$(LIBDIR)

LIBS = $(LIBDIR)/arralloc.o $(LIBDIR)/write_ppm.o $(LIBDIR)/read_options.o

all: $(EXE)

$(OBJ):	$(MAKEFILE) $(SRC) $(INC) utils

.c.o:
	$(MPICC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

utils:
	$(MAKE) -C $(LIBDIR)

$(EXE):	$(OBJ)
	$(MPICC) $(CFLAGS) -o $@ $(OBJ) $(LIBS)


clean: clean-utils
	rm -f $(EXE) *~ $(OBJ) $(EXE) core

clean-utils:
	$(MAKE) -C $(LIBDIR) clean

.PHONY: utils clean-utils
