#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <error.h>
#include <mpi.h>
#include "arralloc.h"
#include "write_ppm.h"
#include "read_options.h"

typedef int (*in_set_fn_t)(const float, const float, const int);

#define EMPTY_TASK -1
#define MIN(X, Y) ((X) < (Y) ? (X) : (Y))

static inline int point_in_mandelbrot_set(const float x0,
                                          const float y0,
                                          const int max_iter)
{
    int i;
    float x;
    float y;
    float x2;
    float y2;
    x = x0;
    y = y0;
    for ( i = 0; i < max_iter; i++ ) {
        x2 = x * x;
        y2 = y * y;
        /* z = (z*z) + c */
        if ( x2 + y2 > 4.0 ) {
            return i;
        } else {
            y = y0 + (2.0 * x * y);
            x = x0 + x2 - y2;
        }
    }
    return max_iter;
}

static inline int point_in_julia_set(const float x0,
                                     const float y0,
                                     const int max_iter)
{
    int i;
    float x;
    float y;
    const float cim = 0.01;
    const float cre = 0.285;
    x = x0;
    y = y0;
    for ( i = 0; i < max_iter; i++ ) {
        if ( x * x + y * y > 4.0 ) {
            return i;
        } else {
            float tmp = y * y;
            y = 2 * x * y + cim;
            x = x * x - tmp + cre;
        }
    }
    return max_iter;
}

/* Initialise data for image array.  Data is stored in "scanline
 * order", i.e. x dimension varies fastest.  You get an array with
 * shape image[grid_size_y][grid_size_x] from this function. */
void initialise_image(int ***image, const int grid_size_x, const int grid_size_y)
{
    int i;
    int j;
    *image = (int**)arralloc(sizeof(int), 2, grid_size_y, grid_size_x);

    if ( NULL == *image ) {
        error(1, errno, "Unable to allocate memory for image\n");
    }
    /* initalise results array to black */
    for ( i = 0; i < grid_size_y; i++ ) {
        for ( j = 0; j < grid_size_x; j++ ) {
            (*image)[i][j] = -1;
        }
    }
}

static inline void send_image_lines(int **image_lines,
                                    int lines[2],
                                    const int grid_size_x,
                                    MPI_Comm comm)
{
    MPI_Send(lines, 2, MPI_INT, 0, 0, comm);
    MPI_Send(&(image_lines[0][0]), grid_size_x * (lines[1] - lines[0]), MPI_INT,
             0, 0, comm);
}

static inline int recv_image_lines(int **image,
                                   const int grid_size_x,
                                   MPI_Comm comm)
{
    int lines[2];
    MPI_Status status;
    /* Receive information about lines from any sender */
    MPI_Recv(lines, 2, MPI_INT, MPI_ANY_SOURCE, 0, comm, &status);
    /* Then get the lines from the same send (using status.MPI_SOURCE
     * to find the source of the previous message) */
    MPI_Recv(&(image[lines[0]][0]), grid_size_x * (lines[1] - lines[0]),
             MPI_INT, status.MPI_SOURCE, 0, comm, &status);

    /* Return where we got the lines from, necessary to send a new task out */
    return status.MPI_SOURCE;
}

int **compute_lines(in_set_fn_t in_set_fn,
                    int lines[2],
                    const float xmin, const float xmax,
                    const float ymin, const float ymax,
                    const int grid_size_x, const int grid_size_y,
                    const int max_iter)
{
    int i;
    int j;
    int start;
    int end;
    int **image_lines;
    float x0;
    float y0;

    start = lines[0];
    end = lines[1];
    initialise_image(&image_lines, grid_size_x, end - start);

    for ( i = 0; i < end - start; i++ ) {
        for ( j = 0; j < grid_size_x; j++ ) {
            /* calculate coordinates x0,y0 of current pixel */
            x0 = xmin + j * ((xmax - xmin) / grid_size_x);
            y0 = ymin + (i + start) * ((ymax - ymin) / grid_size_y);
            image_lines[i][j] = (*in_set_fn)(x0, y0, max_iter);
        }
    }
    return image_lines;
}

/* A task is a set of contiguous lines in the image.  So can be
 * described by an array of two ints.  The task variable holds the
 * current highest numbered line we've asked to be processed. */
static inline void build_task(int lines[2], int *task, int grid_size_y)
{
    const int lines_in_task = 10;
    if ( *task < grid_size_y ) {
        lines[0] = *task;
        lines[1] = MIN(*task + lines_in_task, grid_size_y);
    } else {
        lines[0] = EMPTY_TASK;
        lines[1] = EMPTY_TASK;
    }
    *task += lines_in_task;
}

void master_loop(int **image,
                 const float xmin,
                 const float xmax,
                 const float ymin,
                 const float ymax,
                 const int grid_size_x,
                 const int grid_size_y,
                 const int max_iter,
                 MPI_Comm comm)
{
    int current_task;
    int size;
    int lines[2];
    int i;
    int outstanding_tasks;

    MPI_Comm_size(comm, &size);

    /* Broadcast constant problem data (extent of image, number of
     * pixels, iteration count) */
    {
        float data[4] = {xmin, xmax, ymin, ymax};
        MPI_Bcast(data, 4, MPI_FLOAT, 0, comm);
    }
    {
        int data[3] = {grid_size_x, grid_size_y, max_iter};
        MPI_Bcast(data, 3, MPI_INT, 0, comm);
    }

    outstanding_tasks = 0;
    current_task = 0;
    /* Farm out first round of tasks (all workers are idle) */
    for ( i = 1; i < size; i++ ) {
        build_task(lines, &current_task, grid_size_y);
        MPI_Send(lines, 2, MPI_INT, i, 0, comm);
        /* Don't have to poll for a response if the task was empty */
        if ( lines[0] != EMPTY_TASK )
            outstanding_tasks++;
    }

    /* Sit waiting for data, receive it, then hand out new task */
    while ( outstanding_tasks ) {
        int worker = recv_image_lines(image, grid_size_x, comm);
        --outstanding_tasks;
        build_task(lines, &current_task, grid_size_y);
        MPI_Send(lines, 2, MPI_INT, worker, 0, comm);
        if ( lines[0] != EMPTY_TASK ) {
            ++outstanding_tasks;
        }
    }
}

void worker_loop(in_set_fn_t in_set_fn, MPI_Comm comm)
{
    float xmin, xmax, ymin, ymax;
    int grid_size_x, grid_size_y, max_iter;
    MPI_Status status;
    int lines[2];

    /* Receive constant problem data */
    {
        float data[4];
        MPI_Bcast(data, 4, MPI_FLOAT, 0, comm);
        xmin = data[0];
        xmax = data[1];
        ymin = data[2];
        ymax = data[3];
    }
    {
        int data[3];
        MPI_Bcast(data, 3, MPI_INT, 0, comm);
        grid_size_x = data[0];
        grid_size_y = data[1];
        max_iter = data[2];
    }

    do {
        /* Receive task description */
        MPI_Recv(lines, 2, MPI_INT, 0, 0, comm, &status);
        if ( lines[0] == EMPTY_TASK) {
            /* Nothing do, so return */
            return;
        } else {
            /* Compute the image lines */
            int **data = compute_lines(in_set_fn, lines, xmin, xmax,
                                       ymin, ymax, grid_size_x,
                                       grid_size_y, max_iter);
            send_image_lines(data, lines, grid_size_x, comm);
            free(data);
        }
    } while ( 1 )
        ;
}

void compute_set(in_set_fn_t in_set_fn,
                 int **image,
                 const float xmin,
                 const float xmax,
                 const float ymin,
                 const float ymax,
                 const int grid_size_x,
                 const int grid_size_y,
                 const int max_iter,
                 MPI_Comm comm)
{
    int rank;
    int size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if ( rank == 0 ) {
        master_loop(image, xmin, xmax, ymin, ymax,
                    grid_size_x, grid_size_y, max_iter,
                    comm);
    } else {
        worker_loop(in_set_fn, comm);
    }
}

int main(int argc, char** argv)
{
    int grid_size_x;
    int grid_size_y;
    int max_iter;
    float xmin;
    float xmax;
    float ymin;
    float ymax;
    int **image;
    int rank;
    MPI_Comm comm;
    in_set_fn_t fp;

    MPI_Init(&argc, &argv);

    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);

    read_options(argc, argv, &grid_size_x, &grid_size_y, &max_iter,
                 &xmin, &xmax, &ymin, &ymax);

    fp = &point_in_julia_set;

    if ( rank == 0 ) {
        initialise_image(&image, grid_size_x, grid_size_y);
    }

    compute_set(fp, image,
                xmin, xmax, ymin, ymax,
                grid_size_x, grid_size_y, max_iter,
                comm);

    if ( rank == 0 ) {
        write_ppm("output.ppm", image, grid_size_x, grid_size_y, max_iter);
        free(image);
    }

    MPI_Finalize();
    return 0;
}
