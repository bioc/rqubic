/*
 * Discretize a matrix of expression values (genes in rows, conditions/samples in columns)
 * Author: Jitao David Zhang <jitao_david.zhang@roche.com>
 *
 * C source file to be called by R codes
 */


#include "rqubic.h"

/* quantile function modified from the GNU Scientific Library (GSL) */
static continuous quantile_from_sorted_data (const continuous sorted_data[],
					     const size_t n,
					     const double f) {
  const continuous index = f * (n-1);
  const size_t lhs = (int)index;
  const continuous delta = index - lhs;

  return (1-delta)*sorted_data[lhs]+delta*sorted_data[lhs+1];
}

/* Comparison function for GNU qsort */
static int compare_continuous (const void *a, const void *b)
{
    const continuous *da = a;
    const continuous *db = b;
    
    return (*da > *db) - (*da < *db);
}					 

/* recursive discretization of outliers by quantile
 * both ups and downs are quantiled to discretize outliers
 * the outliers are labelled with -1/1 (the most extreme values) to -rank/rank (the least extreme values)
 * NOTE that the label is reversed as one would assume,i.e. (-rank,-rank+1,...,-1, 0, 1, 2, ..., rank)
 * It is rather (0, -1, -2, ..., -rank, rank, rank-1, ..., 2,1)
 */
discrete discretize_outlier(double value, int rank, 
		    continuous *downs,
		    int cntl,
		    continuous *ups,
		    int cntu) {
  int i;
  double dspace = 1.0/rank;
  for(i=0; i<rank; i++) {
    if (cntl>0 && value<=quantile_from_sorted_data(downs, cntl, dspace * (i+1))) return -i-1;
    if (cntu>0 && value>=quantile_from_sorted_data(ups, cntu,1.0-dspace*(i+1))) return i+1;
  }
  return 0;
}

/* qualitative discretization of the expression matrix */
SEXP discretize_matrix(SEXP exp, SEXP q, SEXP rank){

  int i, j, k, nrow, ncol;
  double *rexp;
  int *rdexp;
  double fhigh, flow, fmed, upper, lower;
  double rq = REAL(q)[0];
  int rrank=INTEGER(rank)[0];
  int cntu, cntl;
  SEXP dim, dexp;
  
  /* serialize matrix */
  rexp=REAL(exp);

  /*get matrix dim*/
  PROTECT(dim = getAttrib(exp, R_DimSymbol));
  nrow = INTEGER(dim)[0];
  ncol = INTEGER(dim)[1];
  PROTECT(dexp = allocMatrix(INTSXP, nrow, ncol));
  rdexp=INTEGER(dexp);

  /* hold row data, no need to initialize, will be filled */
  continuous rowdata[ncol];
  /* hold up/down DEGs, initialize to 0 */
  continuous ups[ncol], downs[ncol];
  memset(ups, 0, ncol*sizeof(*ups));
  memset(downs, 0, ncol*sizeof(*downs));
  
  for(i=0; i<nrow; i++) {
    
    /* fill row data */
    for(j=0; j<ncol; j++)
      rowdata[j] = rexp[i+j*nrow];
    
    qsort(rowdata, ncol, sizeof *rowdata, compare_continuous);
	
    /* define lower and upper boundary of DEG*/
    fhigh = quantile_from_sorted_data(rowdata, ncol, (1 - rq));
    flow = quantile_from_sorted_data(rowdata, ncol, rq);
    fmed = quantile_from_sorted_data(rowdata, ncol, 0.5);
    if((fhigh-fmed) >= (fmed-flow)) {
      upper=2*fmed-flow; lower=flow;
    } else {
      upper=fhigh; lower=2*fmed-fhigh;
    }


    /* fill values of DEGs */
    cntu=0; cntl=0;
    for(k=0; k<ncol; k++) {
      if(rowdata[k] < lower) {downs[cntl]=rowdata[k]; cntl++;}
      if(rowdata[k] > upper) {ups[cntu]=rowdata[k]; cntu++;}
    }
    
    /* discretize outliers */
    for(j=0; j<ncol; j++) 
      rdexp[i+j*nrow]=discretize_outlier(rexp[i+j*nrow],
					 rrank,
					 downs, cntl,
					 ups, cntu);
  }

  /* matrix dim */
  setAttrib(dexp, R_DimNamesSymbol, 
	    getAttrib(exp, R_DimNamesSymbol));

  /* R memory unprotect*/
  UNPROTECT(2);
  return(dexp);
  
}

