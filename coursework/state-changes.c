#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "state-changes.h"
#include "messages.h"
#include "data-structures.h"
#include "helper-functions.h"
#include "parameters.h"

int reproducep(frog f, long * rng_state)
{
    float avg = (float)(f->population_influx) / REPRODUCTION_HOP_COUNT;
    if (f->hop && (f->hop)%REPRODUCTION_HOP_COUNT == 0) {
        return willGiveBirth(avg, rng_state);
    }
    return 0;
}

void maybe_catch_disease(frog f, long * rng_state)
{
    float avg;
    int i;
    /* already infected */
    if (f->diseased) {
        return;
    }

    avg = 0.0f;
    for ( i = 0; i < INFECT_HOP_COUNT; i++ ) {
        avg += (float)(f->infection_level[i]);
    }

    avg /= f->hop < INFECT_HOP_COUNT ? f->hop : INFECT_HOP_COUNT;

    f->diseased = willCatchDisease(avg, rng_state);
}

int diep(frog f, long * rng_state)
{
    return f->hop && f->hop%DIE_HOP_COUNT == 0
        && f->diseased && willDie(rng_state);
}

void update_frog(frog f, int pop_influx, int inf_level)
{
    f->population_influx += pop_influx;
    (f->infection_level)[(f->hop)%INFECT_HOP_COUNT] = inf_level;
    f->hop++;
}

void update_frog_list(list  * frogs, long * rng_state, MPI_Comm comm)
{
    list it;
    int c;
    int buf[3];

    it = *frogs;

    /* Has a control message arrived? In which case we should return
     * immediately and deal with it at the top-level message loop */
    if ( probe_for_interrupt(comm) ) {
        return;
    }
    while ( it ) {
        frog f = (frog)(it->data);
        /* This is not an atomic transaction, the frog can hop cells,
         * but receive an interrupt before sending the message to the
         * relevant cell */
        frogHop(f->x, f->y, &(f->x), &(f->y), rng_state);
        c = getCellFromPosition(f->x, f->y);
        buf[0] = FROG_ARRIVE;
        buf[1] = f->diseased;
        buf[2] = c;
        /* Post the receive that's going to match the send we're about
         * to get from sending a message to the cells.  This can't
         * block for two reasons.  One, we haven't told the cells to
         * send a message yet.  Two, even if we had, the cells might
         * get an interrupt before sending their message.  So we post
         * the receive here, send a message off and clean up later. */
        MPI_Irecv(f->buf, 2, MPI_INT, get_cell_rank(c),
                  FROG_TAG, comm, &(f->req));

        MPI_Send(buf, 3, MPI_INT, get_cell_rank(c),
                 FROG_TAG, comm);

        it = it->next;
    }

    /* clean up outstanding receive requests, but jump to top-level if
     * we receive an interrupt message */
    int not_done;
    do {
        it = *frogs;
        not_done = 0;
        /* If we got a control message, abort to top-level */
        if ( probe_for_interrupt(comm) ) {
            return;
        }
        while ( it ) {
            frog f = (frog)(it->data);
            list tmp = it->next;
            MPI_Status status;
            int flag;
            if ( f->req != NULL ) {
                /* We have an in flight message */
                MPI_Test(&f->req, &flag, &status);
                if ( flag ) {
                    /* It's arrived, update the frog */
                    update_frog(f, f->buf[0], f->buf[1]);
                    if (reproducep(f, rng_state)) {
                        frog new = make_frog(f->x, f->y);
                        push(new, frogs);
                    }
                    maybe_catch_disease(f, rng_state);
                    if (diep(f, rng_state)) {
                        delete_node(frogs, it, &delete_frog);
                    }
                    /* Message is done so reset request status */
                    f->req = NULL;
                } else {
                    /* In flight message has not arrived, so we're not
                     * done here yet */
                    not_done = 1;
                }
            }
            /* This is in case we deleted a frog from the list, in
             * which case it points to the wrong place. */
            it = tmp;
        }
    } while ( not_done )
        ;

}

void update_cell(cell c, int infected)
{
    ++c->population_influx;
    c->infection_level += infected;
}

void update_cell_list(cell_list cells, MPI_Comm comm)
{
    /* Do we have a control message? If so, deal with it at top-level
     * cell loop. */
    if ( probe_for_interrupt(comm) ) {
        return;
    }
    if ( cells->req ) {
        int flag;
        MPI_Status status;
        MPI_Test(&cells->req, &flag, &status);
        if ( flag ) {
            list it;
            cell c;
            for ( it = cells->cells; it; it = it->next ) {
                c = (cell)(it->data);
                if ( cells->buf[0] == FROG_ARRIVE && c->cell == cells->buf[2] ) {
                    /* Does this cell match the one the frog hopped
                     * into? */
                    update_cell(c, cells->buf[1]);
                    cells->buf[0] = c->population_influx - 1;
                    cells->buf[1] = c->infection_level - cells->buf[1];
                    /* This send will always be matched by the
                     * frog (which has just posted a non-block
                     * receive by the time we get here) */
                    MPI_Send(cells->buf, 2, MPI_INT,
                             status.MPI_SOURCE, FROG_TAG,
                             comm);
                }
            }
            MPI_Irecv(cells->buf, 3, MPI_INT, MPI_ANY_SOURCE,
                      FROG_TAG, comm, &cells->req);
        }
    } else {
        MPI_Irecv(cells->buf, 3, MPI_INT, MPI_ANY_SOURCE,
                  FROG_TAG, comm, &cells->req);
    }
}

void reset_cell(cell c)
{
    c->population_influx = 0;
    c->infection_level = 0;
}

void reset_cell_list(cell_list cells)
{
    list it;
    for ( it = cells->cells; it; it = it->next ) {
        reset_cell((cell)(it->data));
    }
}
