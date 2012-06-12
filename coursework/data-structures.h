#ifndef _DATA_STRUCTURES_H
#define _DATA_STRUCTURES_H
#include <mpi.h>
#define REPRODUCTION_HOP_COUNT 300
#define INFECT_HOP_COUNT 500
#define DIE_HOP_COUNT 700

typedef struct {
    float x;
    float y;
    int hop;
    int diseased;
    int population_influx;
    int * infection_level;
    MPI_Request req;
    int buf[3];
} frog_t;

typedef frog_t * frog;

typedef struct {
    int cell;
    int population_influx;
    int infection_level;
} cell_t;

typedef cell_t * cell;

typedef struct list_t {
    void * data;
    struct list_t * prev;
    struct list_t * next;
} list_t;

typedef list_t * list;

typedef struct {
    list cells;
    int buf[3];
    MPI_Request req;
} cell_list_t;

typedef cell_list_t * cell_list;

typedef void (*free_obj_fn)(void *);

frog make_frog(float, float);

cell make_cell(int);

list make_list(void *);

cell_list make_cell_list();

void delete_frog(void *);

void delete_cell(void *);

void delete_node(list *, list, free_obj_fn);

void delete_list(list *, free_obj_fn);

void delete_cell_list(cell_list);

void push(void *, list *);

void print_list(list);

void print_frog(frog);

void print_cell(cell);
#endif
