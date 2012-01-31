#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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
    MPI_Send(&start, 1, MPI_INT, dest, 0, comm);
    MPI_Send(&end, 1, MPI_INT, dest, 0, comm);
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
    int start;
    int end;
    int rank;
    MPI_Comm comm;

    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);

    /* We send these values to master process for writing because we
     * don't assume that the mapping from rank to slice is
     * one-to-one.  For example, we might divide the domain into
     * twice as many slices as there are processes. */
    calc_slice_bounds(slice, nslice, grid_size_y, &start, &end);

    if ( rank == 0 ) {
        int i;
        int j;
        int size;

        /* Copy local bit into global image */
        for ( i = start; i < end; i++ ) {
            for ( j = 0; j < grid_size_x; j++ ) {
                image[i][j] = image_slice[i - start][j];
            }
        }

        /* Receive remote slices */
        MPI_Comm_size(comm, &size);
        for ( i = 1; i < size; i++ ) {
            recv_image_slice(image, grid_size_x, i, comm);
        }
    } else {
        send_image_slice(image_slice, start, end, grid_size_x, 0, comm);
    }
}

int **compute_slice(in_set_fn_t in_set_fn,
                    const int slice, const int nslice,
                    const float xmin, const float xmax,
                    const float ymin, const float ymax,
                    const int grid_size_x, const int grid_size_y,
                    const int max_iter)
{
    int i;
    int j;
    int start;
    int end;
    int **image_slice;
    float x0;
    float y0;

    calc_slice_bounds(slice, nslice, grid_size_y, &start, &end);

    initialise_image(&image_slice, grid_size_x, end - start);

    for ( i = 0; i < end - start; i++ ) {
        for ( j = 0; j < grid_size_x; j++ ) {
            /* calculate coordinates x0,y0 of current pixel */
            x0 = xmin + j * ((xmax - xmin) / grid_size_x);
            y0 = ymin + (i + start) * ((ymax - ymin) / grid_size_y);
            image_slice[i][j] = (*in_set_fn)(x0, y0, max_iter);
        }
    }
    return image_slice;
}

void compute_set(in_set_fn_t in_set_fn,
                 int **image,
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
    image_slice = compute_slice(in_set_fn, slice, nslice, xmin, xmax, ymin, ymax,
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
    in_set_fn_t fp;

    MPI_Init(&argc, &argv);

    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);

    read_options(argc, argv, &grid_size_x, &grid_size_y, &max_iter,
                 &xmin, &xmax, &ymin, &ymax);

    if ( rank == 0 ) {
        initialise_image(&image, grid_size_x, grid_size_y);
    }

    fp = &point_in_julia_set;

    compute_set(fp, image,
                xmin, xmax, ymin, ymax,
                grid_size_x, grid_size_y, max_iter);

    if ( rank == 0 ) {
        write_ppm("output.ppm", image, grid_size_x, grid_size_y, max_iter);
        free(image);
    }

    MPI_Finalize();
    return 0;
}
