#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <error.h>
#include <mpi.h>
#include "arralloc.h"
#include "write_ppm.h"
#include "read_options.h"

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

/* Divide up Y dimension to calculate bounds of this slice */
static inline void calc_slice_bounds(const int slice, const int nslice,
                                     const int grid_size_y,
                                     int *start, int *end)
{
    *start = slice * grid_size_y / nslice;
    *end = (slice + 1) * grid_size_y / nslice;
}

static inline void send_image_slice(int **image_slice,
                                    int start,
                                    int end,
                                    const int grid_size_x,
                                    int dest,
                                    MPI_Comm comm)
{
    MPI_Send(&start, 1, MPI_INT, 0, 0, comm);
    MPI_Send(&end, 1, MPI_INT, 0, 0, comm);
    MPI_Send(&(image_slice[0][0]), grid_size_x * (end - start), MPI_INT,
             dest, 0, comm);
}

static inline void recv_image_slice(int **image,
                                    const int grid_size_x,
                                    int source,
                                    MPI_Comm comm)
{
    int start;
    int end;
    MPI_Status status;
    MPI_Recv(&start, 1, MPI_INT, source, 0, comm, &status);
    MPI_Recv(&end, 1, MPI_INT, source, 0, comm, &status);
    MPI_Recv(&(image[start][0]), grid_size_x * (end - start), MPI_INT,
             source, 0, comm, &status);
}

void copy_slice_to_image(int **image_slice, int **image,
                         const int slice, const int nslice,
                         const int grid_size_x, const int grid_size_y)
{
    int rank;
    MPI_Comm comm;

    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);

    if ( rank == 0 ) {
        /* Copy local slices into global image */
        /* Receive remote slices into global image */
    } else {
        /* Send slice to rank 0 */
    }
}

int **compute_mandelbrot_slice(const int slice, const int nslice,
                               const float xmin, const float xmax,
                               const float ymin, const float ymax,
                               const int grid_size_x, const int grid_size_y,
                               const int max_iter)
{
    /*
     * Compute the mandelbrot set in a SLICE of the whole domain.
     * The total number of slices is given by NSLICE.
     *
     * The utility function calc_slice_bounds is provided to work out
     * the bounds of a slice.  It divides the image in the y
     * direction.
     *
     * This function should return a newly initialised image slice.
     */
    return NULL;
}

void compute_mandelbrot_set(int **image,
                            const float xmin,
                            const float xmax,
                            const float ymin,
                            const float ymax,
                            const int grid_size_x,
                            const int grid_size_y,
                            const int max_iter)
{
    int rank;
    int size;
    int slice;
    int nslice;
    int **image_slice;
    MPI_Comm comm;

    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    nslice = size;
    slice = rank;
    /* Compute local image slice */
    image_slice = compute_mandelbrot_slice(slice, nslice, xmin, xmax, ymin, ymax,
                                           grid_size_x, grid_size_y, max_iter);

    /* Copy into global image, note IMAGE is only defined on rank 0 */
    copy_slice_to_image(image_slice, image, slice, nslice, grid_size_x, grid_size_y);

    free(image_slice);
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

    MPI_Init(&argc, &argv);

    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);

    read_options(argc, argv, &grid_size_x, &grid_size_y, &max_iter,
                 &xmin, &xmax, &ymin, &ymax);

    if ( rank == 0 ) {
        initialise_image(&image, grid_size_x, grid_size_y);
    }

    compute_mandelbrot_set(image, xmin, xmax, ymin, ymax,
                           grid_size_x, grid_size_y, max_iter);

    if ( rank == 0 ) {
        write_ppm("output.ppm", image, grid_size_x, grid_size_y, max_iter);
        free(image);
    }

    MPI_Finalize();
    return 0;
}
