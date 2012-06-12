CC = mpicc
CFLAGS = -Wall -Wextra -O3
LDFLAGS =

OBJ = data-structures.o helper-functions.o ran2.o messages.o state-changes.o

EXE = main

all: $(EXE)

main: main.c $(OBJ)

data-structures.o: data-structures.c data-structures.h helper-functions.h

helper-functions.o: helper-functions.c helper-functions.h ran2.h

ran2.o: ran2.c ran2.h

messages.o: messages.c messages.h data-structures.h state-changes.h helper-functions.h parameters.h

state-changes.o: state-changes.c data-structures.h helper-functions.h parameters.h messages.h

clean:
	-rm -f $(EXE) $(OBJ)

.PHONY: clean
