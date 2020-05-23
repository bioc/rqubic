#ifndef _RQUBIC_H
#define _RQUBIC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <ctype.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Compatibility of __attribute__ with non-GNU */
#ifndef __GNUC__
#  define __attribute__(x) /* Nothing */
#endif

/* Pretend that C has boolean type */
#define TRUE 1
#define FALSE 0
#define boolean unsigned char
#ifndef __cplusplus
#ifndef bool
#define bool unsigned char
#endif
#endif

/* Macros */
#define MAX(a,b)  ((a)>(b)?(a):(b))
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define ABS(x)    ((x)>0?(x):-(x))

/* Variable and array allocation */
#define AllocVar(pt) (pt = xmalloc(sizeof(*pt)))
#define AllocArray(pt, size) (pt = xmalloc(sizeof(*pt) * (size)))
#define ReAllocArray(pt, size) (pt = xrealloc(pt, sizeof(*pt) * (size)))

/* Two major data types */
typedef double continuous;
typedef int discrete;
typedef unsigned short int bits16;

/* Wrapper for memory allocations */
void *xmalloc ( int size );
/* Wrapper for memory re-allocations */
void *xrealloc ( void* ptr, int size );

/* structs used by more than one module */
typedef struct Edge{
	int gene_one;
	int gene_two;
	int score;
} Edge, *pEdge, **ppEdge;

/* initialization */
// names
extern SEXP RQUBIC_edgelist_tag;

SEXP RQUBIC_init(SEXP env);

static const char *EDGELIST_ATT_NAME="edgelist";
static const char *MINCOL_ATT_NAME="minimumCol";


/* public functions*/
SEXP discretize_matrix(SEXP exp, SEXP q, SEXP rank);
SEXP generate_sorted_seeds(SEXP dexp, SEXP col_width);
SEXP qubicluster(SEXP seeds, SEXP exprs, SEXP sigma_val, SEXP symbols_val, SEXP report_no,
			SEXP tolerance_val, SEXP filter_proportion);
#endif
