#ifndef _BICLUSTER_H
#define _BICLUSTER_H

#include "rqubic.h"

/* structs */
/* dynamically allocated stack */
typedef struct dyStack{
  int top;             /* top element index */
  int items[];		   /* data storage */
} dyStack;


/* this structure holds the matching score for each row */
struct rowMatch{
    int row_id;
    int matches;
};

typedef struct Block{
	struct dyStack *genes;
	struct dyStack *conds;
	int score;
	int block_rows;
	int block_cols;
	int block_rows_pre;
	double significance;
} Block;

/* prototypes */

// seed related
void seed_intersect (const discrete *s1, const discrete *s2, 
		     const int ncol, bits16 **profile);
void seed_update (const discrete *s, const int ncol, bits16 **profile);
void seed_deduct (const discrete *s, const int ncol, discrete **profile);
void seed_current (const discrete *s, int* cnt, int components,
		   const int ncol, const int sigma, bits16 **profile);
void seed_intersect_r (const discrete *s, int *cnt, int *m_cnt, int components,
		       const int ncol, const int sigma, const double tolerance,
		       const bits16 **profile);


SEXP cluster(SEXP seeds, SEXP exprs, SEXP sigma_val, SEXP symbols_val,
	     SEXP report_no, SEXP tolerance, SEXP filter_proportion);
static void sort_block_list(Block **el, int n);
static int intersect_row(const bool *colcand, const discrete *g1, const discrete *g2, 
			 const int ncol);
static int reverse_row(const bool *colcand, const discrete *g1, const discrete *g2,
		       const int ncol, const int *symbols);
void scan_block (struct dyStack *gene_set, Block *b_ptr,
		 discrete **arr_c,const int ncol,
		 const int sigma, const double tolerance, bits16 **profile);

/* Stack-related operations */
bool isInStack(struct dyStack *ds, int element);
bool isBothInStack(struct dyStack *ds, int element1, int element2);

struct dyStack *dsNew(const int size);
void dsPush(struct dyStack *ds, int element);
int dsIntersect(struct dyStack *ds1, struct dyStack *ds2);
#define dsItem(pds, j) ((pds)->items[j])
#define dsSize(pds) ((pds)->top + 1)
#define dsClear(pds) ((pds)->top = -1)


SEXP report_blocks(Block** bb, 
		  const int num, 
		  const int reportBlockNo,
		  const double filter_proportion,
		  discrete **arr_c, const int *symbols);
/* profile: 2-dimension unsigned short int array
 * for most C implementations, bits16 should be typedef'ed to
 * unsigned short int. It has been explicitly defined in the rqubic.h
*/
#endif
