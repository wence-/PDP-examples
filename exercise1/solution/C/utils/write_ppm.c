#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <error.h>
#include <math.h>
#include "write_ppm.h"

/*
 * Convert a single iteration count to an RGB value by treating it as
 * an HSV triplet with saturation and value both set as 1 
 */
static inline void gray_to_rgb(int gray, int rgb[3])
{
    double h;
    double s = 1;
    double v = MAX_COLOUR_VALS;
    double f, p, q, t;

    h = (360.0 * gray) / (60 * MAX_COLOUR_VALS);
    if ( h < 0 ) {
        /* Invalid colour, set to black */
        rgb[R] = 0;
        rgb[G] = 0;
        rgb[B] = 0;
        return;
    }
    f = h - (int)h;
    
    p = v * (1 - s);
    q = v * (1 - s * f);
    t = v * (1 - s * (1 - f));
    switch ( (int)h ) {
    case 0:
        rgb[R] = (int)v;
        rgb[G] = (int)t;
        rgb[B] = (int)p;
        break;
    case 1:
        rgb[R] = (int)q;
        rgb[G] = (int)v;
        rgb[B] = (int)p;
        break;
    case 2:
        rgb[R] = (int)p;
        rgb[G] = (int)v;
        rgb[B] = (int)t;
        break;
    case 3:
        rgb[R] = (int)p;
        rgb[G] = (int)q;
        rgb[B] = (int)v;
        break;
    case 4:
        rgb[R] = (int)t;
        rgb[G] = (int)p;
        rgb[B] = (int)v;
        break;
    case 5:
        rgb[R] = (int)v;
        rgb[G] = (int)p;
        rgb[B] = (int)q;
        break;
    default:
        rgb[R] = (int)0;
        rgb[G] = (int)0;
        rgb[B] = (int)0;
    }
}

/*
 * Write a PPM file containing the NI x NJ pixels of IMAGE to FILE.
 *
 * IMAGE is considered as a set of grayscale values that are converted
 * to RGB by mapping them onto HSV.
 */
void write_ppm(const char *file, int **image, int ni, int nj)
{
    int i;
    int j;
    FILE *fp;
    int rgb[3];
    fp = fopen(file,"w");

    if ( NULL == fp ) {
        error(1, errno, "Unable to open file %s for output\n", file);
    }

    /*
     * PPM format is:
     * P3
     * WIDTH HEIGHT
     * MAX_COLOURS
     * R G B
     * R G B
     * R G B
     * ...
     */
    fprintf(fp, "P3\n%d %d\n%d\n", ni, nj, MAX_COLOUR_VALS);

    for ( i = 0; i < ni; i++ ) {
        for ( j = 0; j < nj; j++ ) {
            gray_to_rgb(image[i][j], rgb);
            fprintf (fp,"%d %d %d\n", rgb[R], rgb[G], rgb[B]);
        }
    }

    fclose(fp);
}
