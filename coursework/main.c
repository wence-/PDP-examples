#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "data-structures.h"
#include "messages.h"
#include "parameters.h"
#include "helper-functions.h"

int main(int argc, char **argv)
{
    MPI_Comm comm;
    int rank;
    int size;
    long rng_state;
    MPI_Init(&argc, &argv);

    comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if ( size < 4 ) {
        if ( !rank ) {
            fprintf(stderr, "Unable to run with fewer than 4 MPI ranks\n");
        }
        MPI_Finalize();
        return -1;
    }

    rng_state = -1 - rank;
    initialiseRNG(&rng_state);
    if ( 0 == rank ) {
        master_spin(size, comm);
    } else if ( rank < 3 ) {
        cell_spin(rank, NUM_CELLS, comm);
    } else {
        frog_spin(rank, size, NUM_FROGS, &rng_state, comm);
    }

    MPI_Finalize();
    return 0;
}
