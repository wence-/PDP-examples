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

static inline int point_in_mandelbrot_set(const float x0,
                                          const float y0,
                                          const int max_iter)
{
    return -1;
}

static inline int point_in_julia_set(const float x0,
                                     const float y0,
                                     const int max_iter)
{
    return -1;
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
}

void worker_loop(in_set_fn_t in_set_fn, MPI_Comm comm)
{
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
