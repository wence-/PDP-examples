#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <error.h>
#include "arralloc.h"
#include "write_ppm.h"
#include "read_options.h"

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

void copy_slice_to_image(int **image_slice, int **image,
                         const int slice, const int nslice,
                         const int grid_size_x, const int grid_size_y)
{
    int i;
    int j;
    int start;
    int end;

    calc_slice_bounds(slice, nslice, grid_size_y, &start, &end);

    for ( i = start; i < end; i++ ) {
        for ( j = 0; j < grid_size_x; j++ ) {
            image[i][j] = image_slice[i - start][j];
        }
    }
}

int **compute_mandelbrot_slice(const int slice, const int nslice,
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
            image_slice[i][j] = point_in_mandelbrot_set(x0, y0, max_iter);
        }
    }
    return image_slice;
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
    int slice;
    int nslice;
    int **image_slice;

    /* Arbitrary number of slices */
    nslice = 3;
    for ( slice = 0; slice < nslice; slice++ ) {
        image_slice = compute_mandelbrot_slice(slice, nslice, xmin, xmax,
                                               ymin, ymax,
                                               grid_size_x, grid_size_y,
                                               max_iter);
        copy_slice_to_image(image_slice, image, slice, nslice,
                            grid_size_x, grid_size_y);
        free(image_slice);
    }
}

int main(int argc, char** argv)
{
    int grid_size_x = GRIDSIZE_X;
    int grid_size_y = GRIDSIZE_Y;
    int max_iter = ITERATIONS;
    float xmin = XMIN;
    float xmax = XMAX;
    float ymin = YMIN;
    float ymax = YMAX;
    int **image;

    read_options(argc, argv, &grid_size_x, &grid_size_y, &max_iter,
                 &xmin, &xmax, &ymin, &ymax);

    initialise_image(&image, grid_size_x, grid_size_y);

    compute_mandelbrot_set(image, xmin, xmax, ymin, ymax,
                           grid_size_x, grid_size_y, max_iter);

    write_ppm("output.ppm", image, grid_size_x, grid_size_y, max_iter);

    free(image);
    return 0;
}
