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
        if ( x2 + y2 >= 4.0 ) {
            return i;
        } else {
            y = y0 + (2.0 * x * y);
            x = x0 + x2 - y2;
        }
    }
    return max_iter;
}

void compute_mandelbrot_set(int **image,
                            const float xmin,
                            const float xmax,
                            const float ymin,
                            const float ymax,
                            const int grid_size,
                            const int max_iter)
{
    int i;
    int j;
    float x0;
    float y0;
    for ( i = 0; i < grid_size; i++ ) {
        for ( j = 0; j < grid_size; j++ ) {
            /* calculate coordinates x0,y0 of current pixel */
            x0 = xmin + i * ((xmax - xmin) / grid_size);
            y0 = ymin + j * ((ymax - ymin) / grid_size);
            image[i][j] = point_in_mandelbrot_set(x0, y0, max_iter);
        }
    }

}

void initialise_image(int ***image, const int grid_size)
{
    int i;
    int j;
    *image = (int**)arralloc(sizeof(int), 2, grid_size, grid_size);

    if ( NULL == *image ) {
        error(1, errno, "Unable to allocate memory for image\n");
    }
    /* initalise results array to black */
    for ( i = 0; i < grid_size; i++ ) {
        for ( j = 0; j < grid_size; j++ ) {
            (*image)[i][j] = -1;
        }
    }
}

int main(int argc, char** argv)
{
    int grid_size = GRIDSIZE;
    int max_iter = ITERATIONS;
    float xmin = XMIN;
    float xmax = XMAX;
    float ymin = YMIN;
    float ymax = YMAX;
    int **image;

    read_options(argc, argv, &grid_size, &max_iter, &xmin, &xmax, &ymin, &ymax);

    initialise_image(&image, grid_size);
    write_ppm("input.ppm", image, grid_size, grid_size);

    compute_mandelbrot_set(image, xmin, xmax, ymin, ymax, grid_size, max_iter);

    /* Sort out colours and write out the output file. */
    write_ppm("output.ppm", image, grid_size, grid_size);

    free(image);
    return 0;
}
