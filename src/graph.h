#ifndef _GRAPH_H
#define _GRAPH_H

#include "rqubic.h"
#include "fib.h"

/* prototypes */
static int str_intersect_r (const discrete *s1, const discrete *s2, const int ncol);

/* static int compare_edges (const void *a, const void *b);*/
static void fh_insert_fixed (struct fibheap *a, Edge *i, Edge **cur_min);

static void edgelistFinalizer (SEXP ptr);

#endif
