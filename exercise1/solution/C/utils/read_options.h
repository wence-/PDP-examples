#ifndef READ_OPTIONS_H
#define READ_OPTIONS_H

void read_options(int, char**, int*, int*, float*, float*, float*, float*);

/* size of entire GRID (in one dimenson)*/
#define GRIDSIZE 768


/* default coordinates for calculation */
#define XMIN (0.32)
#define XMAX (0.39)
#define YMIN (0.37)
#define YMAX (0.44)

/* 
 * #define X 0.0
 * #define Y 0.0
 */

/* default number of iterations */
#define ITERATIONS 5000

#endif
