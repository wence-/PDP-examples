#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#include "messages.h"
#include "data-structures.h"
#include "state-changes.h"
#include "helper-functions.h"
#include "parameters.h"

int probe_for_interrupt(MPI_Comm comm)
{
    int flag;
    MPI_Status status;
    MPI_Iprobe(0, CONTROL_TAG, comm, &flag, &status);
    return flag;
}

void master_spin(int size, MPI_Comm comm)
{
    int i;
    int j;
    int msg;
    i = 0;
    /* Sleep for a specified amount then wake and send an end of year
     * message to everyone.  After a specified number of "years"
     * indicate that the simulation should end. */
    do {
        usleep(500000);
        msg = YEAR_END;
        for ( j = 1; j < size; j++ ) {
            MPI_Send(&msg, 1, MPI_INT, j, CONTROL_TAG, comm);
        }
    } while ( i++ < 10 )
        ;
    msg = SIMULATION_END;
    for ( j = 1; j < size; j++ ) {
        MPI_Send(&msg, 1, MPI_INT, j, CONTROL_TAG, comm);
    }
}

void cell_spin(int rank, int ncell, MPI_Comm comm)
{
    cell_list cells;
    int i;
    int msg;
    cells = make_cell_list();
    /* Build the list of cells */
    for ( i = 0; i < ncell; i++ ) {
        if ( rank == get_cell_rank(i) ) {
            push(make_cell(i), &(cells->cells));
        }
    }
    i = 0;
    /* Loop, updating the list of cells.  Every time through the loop,
     * we look to see if we have a control message (from the master)
     * and if so, act accordingly (reset cells for end of year, or
     * finish simulation).  Note the the control message could arrive
     * while in the update_cell_list function, so that also probes for
     * interrupts and returns early if it gets one. */
    do {
        MPI_Status status;
        int flag;
        msg = -1;
        update_cell_list(cells, comm);
        MPI_Iprobe(0, CONTROL_TAG, comm, &flag, &status);
        if ( flag ) {
            MPI_Recv(&msg, 1, MPI_INT, status.MPI_SOURCE,
                     status.MPI_TAG, comm, &status);
        }
        if ( msg == YEAR_END ) {
            list it;
            for ( it = cells->cells; it; it = it->next ) {
                print_cell((cell)(it->data));
            }
            reset_cell_list(cells);
        }
    } while ( msg != SIMULATION_END )
        ;

    delete_cell_list(cells);
}

void frog_spin(int rank, int size, int nfrog, long * rng_state, MPI_Comm comm)
{
    list frogs = NULL;
    int i;
    int start;
    int end;
    int chunk;
    int msg;

    chunk = nfrog/(size - 3);
    start = chunk * (rank - 3);
    end = start + chunk;
    if ( rank == size - 1 ) {
        end = nfrog;
    }
    /* Build initial list of frogs */
    for ( i = start; i < end; i++ ) {
        frog f = make_frog(0,0);
        frogHop(0, 0, &f->x, &f->y, rng_state); /* random starting position */
        push(f, &frogs);
    }
    i = 0;
    /* Loop as for cell list, looking for interrupts.  To ensure we
     * don't get deadlocks, lots of magic happens in update_frog_list */
    do {
        MPI_Status status;
        int flag;
        msg = -1;
        update_frog_list(&frogs, rng_state, comm);
        MPI_Iprobe(0, CONTROL_TAG, comm, &flag, &status);
        if ( flag ) {
            MPI_Recv(&msg, 1, MPI_INT, status.MPI_SOURCE,
                     status.MPI_TAG, comm, &status);
        }
        if ( msg == YEAR_END ) {
            int nfrogs[2] = {0,0};
            list it = frogs;
            /* Calculate global number of frogs (and diseased
             * frogs) */
            while ( it ) {
                nfrogs[0]++;
                nfrogs[1] += ((frog)(it->data))->diseased;
                it = it->next;
            }
            /* Should do this by building a sub-communicator and doing
             * a reduction, but that's hard work */
            if ( rank == 3 ) {
                int tmp[2];
                int j;
                for ( j = 4; j < size; j++ ) {
                    MPI_Recv(&tmp, 2, MPI_INT, j, FROG_REDUCTION, comm, &status);
                    nfrogs[0] += tmp[0];
                    nfrogs[1] += tmp[1];
                }
                printf("There are %d frogs (%d infected) in year %d\n",
                       nfrogs[0], nfrogs[1], i++);
            } else {
                MPI_Send(&nfrogs, 2, MPI_INT, 3, FROG_REDUCTION, comm);
            }
        }
    } while ( msg != SIMULATION_END )
        ;
    delete_list(&frogs, &delete_frog);
}
