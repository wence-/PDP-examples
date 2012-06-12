#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <mpi.h>
#include "data-structures.h"
#include "helper-functions.h"

frog make_frog(float x, float y)
{
    frog ret = (frog)malloc(sizeof(frog_t));
    ret->x = x;
    ret->y = y;
    ret->diseased = 0;
    ret->population_influx = 0;
    ret->hop = 0;
    ret->infection_level = (int *)calloc(INFECT_HOP_COUNT, sizeof(int));
    ret->req = NULL;
    return ret;
}

cell make_cell(int c)
{
    cell ret = (cell)malloc(sizeof(cell_t));
    ret->cell = c;
    ret->population_influx = 0;
    ret->infection_level = 1;
    return ret;
}

list make_list(void * data)
{
    list ret = (list)malloc(sizeof(list_t));
    ret->data = data;
    ret->prev = NULL;
    ret->next = NULL;
    return ret;
}

cell_list make_cell_list()
{
    cell_list ret = (cell_list)malloc(sizeof(cell_list_t));
    ret->cells = NULL;
    ret->buf[0] = -1;
    ret->buf[1] = -1;
    ret->buf[2] = -1;
    ret->req = NULL;
    return ret;
}

void delete_frog(void * data)
{
    frog x = (frog)data;
    if (x) {
        free(x->infection_level);
        free(x);
    }
    x = NULL;
}

void delete_cell(void * data)
{
    cell x = (cell)data;
    if (x) {
        free(x);
    }
    x = NULL;
}

void delete_node(list * ll, list node, free_obj_fn free_fn)
{
    if ( *ll == NULL || node == NULL ) {
        return;
    }

    if ( *ll == node ) {
        (*ll) = node->next;
    }

    if ( node->next ) {
        node->next->prev = node->prev;
    }

    if ( node->prev ) {
        node->prev->next = node->next;
    }

    (*free_fn)(node->data);
    free(node);
}

void delete_cell_list(cell_list cells)
{
    delete_list(&(cells->cells), &delete_cell);
    free(cells);
}

void delete_list(list * ll, free_obj_fn free_fn)
{
    list l = *ll;
    if (!l) {
        return;
    }
    while (l->next) {
        l = l->next;
        (*free_fn)(l->prev->data);
        free(l->prev);
    }
    (*free_fn)(l->data);
    free(l);
    *ll = NULL;
}

void push(void * data, list * l)
{
    list head = make_list(data);
    if (*l) {
        head->next = *l;
        (*l)->prev = head;
    }
    *l = head;
}

void print_frog(frog f)
{
    if (f) {
        printf("frog at (%d) hop %d (infection %d)\n",
               getCellFromPosition(f->x, f->y),
               f->hop,
               f->diseased);
    }
}

void print_cell(cell c)
{
    if (c) {
        printf("%d %d\n", c->population_influx, c->infection_level);
    }
}

void print_list(list l)
{
    while (l) {
        printf("<%p>: data <%p>; prev <%p>; next <%p>\n",
               l, l->data, l->prev, l->next);
        l = l->next;
    }
}
