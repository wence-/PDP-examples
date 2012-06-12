#ifndef _MESSAGES_H
#define _MESSAGES_H

#include <mpi.h>
#include "data-structures.h"
#include "state-changes.h"

int probe_for_interrupt(MPI_Comm);
void master_spin(int, MPI_Comm);
void cell_spin(int, int, MPI_Comm);
void frog_spin(int, int, int, long *, MPI_Comm);
#endif
