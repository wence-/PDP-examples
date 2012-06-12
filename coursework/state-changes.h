#ifndef _STATE_CHANGES_H
#define _STATE_CHANGES_H
#include "data-structures.h"

int reproducep(frog, long *);

void maybe_catch_disease(frog, long *);

int diep(frog, long *);

void update_frog(frog, int, int);

void update_frog_list(list *, long *, MPI_Comm comm);

void update_cell(cell, int);

void update_cell_list(cell_list, MPI_Comm comm);

void reset_cell(cell);

void reset_cell_list(cell_list);

static inline int get_cell_rank(int c)
{
    return c < 8 ? 1 : 2;
}

#endif
