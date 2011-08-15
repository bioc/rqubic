/*
 * Establish weighted graph model from quantileDiscretized expression matrix
 * Generating seeds and identifying biclusters
 * 
 * C source file to be called by R codes
 */

#include "graph.h"

/* Note: reduce the HEAP_SIZE when the data contains so many genes that memory is not enough*/

static const int HEAP_SIZE = 20000000;

/**************************************************************/
/*calculate the weight of edges in the complete graph */
static int str_intersect_r (const discrete *s1, const discrete *s2, const int ncol) {
  int common_cnt = 0;
  int i;
  for (i=0; i<ncol; i++) {
    if (*s1 == *s2 && (*s1!=0)) common_cnt++;
    s1++;s2++;
  }
  return common_cnt;
}


/**************************************************************/
/* Fibonacci heap related subroutine */
/* compare edges by scores */
static int compare_edges(void *a, void *b) {
  const int score_a = ((pEdge)a)->score;
  const int score_b = ((pEdge)b)->score;

  return (score_a > score_b) - (score_a < score_b);
}

static void fh_insert_fixed(struct fibheap *a, pEdge i, ppEdge cur_min) {
  if (a->fh_n < HEAP_SIZE) {
    fh_insert(a, (void *)i);
  } else {
    if (compare_edges(cur_min, i) < 0) {
      /* Remove least value and renew */
      fh_extractmin(a);
      fh_insert(a, (void *)i);
      /* Keep a memory of the current min */
      *cur_min = (pEdge)fh_min(a);
    }
  }
}

static void fh_dump(struct fibheap *a, ppEdge res) {
    int i;
    int n = a->fh_n;
    for (i=n-1; i>=0; i--)
        res[i] = (pEdge) fh_extractmin(a);
}

static void *sexp2ptr(SEXP s, Rboolean null_ok, SEXP tag, char *type) {
  void *p;
  if (TYPEOF(s) != EXTPTRSXP || R_ExternalPtrTag(s) != tag)
    error("bad %s pointer", type);
  p = R_ExternalPtrAddr(s);
  if (!null_ok && p == NULL)
    error("null %s pointer", type);
  return(p);
}

static void edgelistFinalizer(SEXP ptr) {
  if(!R_ExternalPtrAddr(ptr)) return;
  Edge** edge_list = sexp2ptr(ptr, FALSE, RQUBIC_edgelist_tag, "ppEdge");
  int i=0;
  /* Note that external pointer protected field points to the array length*/
  int rec_num=INTEGER(R_ExternalPtrProtected(ptr))[0];

  /* to get the rec_num, one can also set an attribute */
  /* int rec_num=INTEGER(getAttrib(ptr, install(EDGELIST_ATT_NAME)))[0];*/

  for(i=0; i<rec_num; i++) {
    free(edge_list[i]);
  }
  free(edge_list);
  R_ClearExternalPtr(ptr);
}

/**************************************************************/

/* build weighted graphs */
SEXP generate_sorted_seeds (SEXP dexp, SEXP col_width) {
  int nrow, ncol;
  int i,j, cnt;
  int rec_num=0;
  int mincol = INTEGER(col_width)[0];
  int *ddexp;
  discrete **arr_c;
  ppEdge edge_list;
  pEdge edge_ptr;
  SEXP dim, ptr, class, ans;

  /*get matrix dim*/
  PROTECT(dim = getAttrib(dexp, R_DimSymbol));
  nrow = INTEGER(dim)[0];
  ncol = INTEGER(dim)[1];
  UNPROTECT(1);

  /* serialize matrix into a 2d-array
   * note that in QUBIC, arr_c should be a index of levels, but not the levels themselves
   * there, symbols[arr_c[i][j]] gives the "real level"
   * here, we use arr_c[i][j] to save the real level
  */

  ddexp = INTEGER(dexp);
  AllocArray(arr_c, nrow);
  for(i=0; i<nrow; i++) {
    AllocArray(arr_c[i], ncol);
    for(j=0; j<ncol; j++) {
      arr_c[i][j] = ddexp[i + nrow * j];  /* Attention to the index of matrix */
    }
  }
  /*minimum column width of a bicluster*/
  if(mincol==2) mincol=MAX(ncol/20, 2);

  AllocArray(edge_list, HEAP_SIZE);
  
  /* Allocating heap structure and set compare functions*/
  struct fibheap *heap;
  heap = fh_makeheap();
  fh_setcmp(heap, compare_edges);
  
  /* set current min edge, and references */
  Edge __cur_min = {0, 0, mincol};
  pEdge _cur_min = &__cur_min;
  ppEdge cur_min = &_cur_min;

  /* iterate over all genes to retrieve all edges 
   * and add to the heap if the heap is not yet full or it is larger than the minimum score 
  */
  for(i=0;i<nrow;i++) {
    for(j=i+1;j<nrow;j++) {
      cnt = str_intersect_r(arr_c[i], arr_c[j], ncol);
      if (cnt < (*cur_min)->score) continue;

      AllocVar(edge_ptr);
      edge_ptr->gene_one=i;
      edge_ptr->gene_two=j;
      edge_ptr->score=cnt;

      fh_insert_fixed(heap, edge_ptr, cur_min);
    }
  }

  rec_num = heap->fh_n;
  if(rec_num==0) {
    fprintf(stderr, "No enough overlap between genes\n");
    return R_NilValue;
  }
    
  /* sort the seeds (highest to lowest) */
  ReAllocArray(edge_list, rec_num);
  fh_dump(heap, edge_list);

  /* Make external pointer of the edge list and return */
  /* Note that the prot pointer is used to save the length of dynamic array here*/
  ptr = R_MakeExternalPtr(edge_list, RQUBIC_edgelist_tag, ScalarInteger(rec_num));
  /* to set the rec_num, one can also set an attribute*/
  /* setAttrib(ptr, install(EDGELIST_ATT_NAME), ScalarInteger(rec_num)); */
  PROTECT(ptr);
  R_RegisterCFinalizerEx(ptr, edgelistFinalizer, TRUE);
    
  PROTECT(ans = allocVector(INTSXP, 1));
  INTEGER(ans)[0] = rec_num;
  setAttrib(ans, install(EDGELIST_ATT_NAME), ptr);
  setAttrib(ans, install(MINCOL_ATT_NAME), ScalarInteger(mincol));

  PROTECT(class=allocVector(STRSXP, 1));
  SET_STRING_ELT(class, 0, mkChar("rqubicSeeds"));
  classgets(ans, class);

  UNPROTECT(3);
  return ans;

  /*
   *  blocks_number = cluster(edge_list, rec_num);
   *
  */

}

