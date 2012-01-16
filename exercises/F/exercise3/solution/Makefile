# Leave alone below here


#
# lists of sources, etc
#
SRC = fractal-set.f90

EXE = fractal-set

include ../../../config.mak

LIBDIR = ../../utils
FFLAGS := $(FFLAGS) -I$(LIBDIR)

LIBS = $(LIBDIR)/write_ppm.o $(LIBDIR)/read_options.o

all: $(EXE)

utils:
	$(MAKE) -C $(LIBDIR)

$(EXE):	$(SRC) utils
	$(MPIF90) $(FFLAGS) -o $@ $(SRC) $(LIBS)


clean: clean-utils
	rm -f $(EXE) *~ $(OBJ) $(EXE) core

clean-utils:
	$(MAKE) -C $(LIBDIR) clean

.PHONY: utils clean-utils
