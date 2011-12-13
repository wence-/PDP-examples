#ifndef READ_OPTIONS_H
#define READ_OPTIONS_H

void read_options(int, char**, int*, int*, int *,
                  float*, float*, float*, float*);


/* default coordinates for calculation */
#define XMIN (-2.0)
#define XMAX (1.0)
#define YMIN (-1.2)
#define YMAX (1.2)

/* size of entire GRID */
#define GRIDSIZE_X 768
#define GRIDSIZE_Y (int)(GRIDSIZE_X * (YMAX - YMIN)/(XMAX - XMIN))

/* default number of iterations */
#define ITERATIONS 5000

#endif
