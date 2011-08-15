/* Bicluster algorithm of QUBIC
 */

#include "bicluster.h"

/**************************************************************/
static void block_init(Edge *e, Block *b, 
		       struct dyStack *genes, struct dyStack *scores,
		       bool *candidates, const int cand_threshold,
		       int *components, struct dyStack *allincluster,
		       const int nrow, const int ncol, discrete ** const arr_c, const int mincol);


/* block */
static int compare_int (const void *a, const void *b) {
  const int *da = a;
  const int *db = b;
  return (*da > *db)-(*da < *db);
}

static void update_colcand(bool *colcand, const discrete *g1, const discrete *g2,
			   const int ncol) {
  int i;
  for(i=0; i<ncol;i++) 
    if(colcand[i] && (g1[i] != g2[i]))
      colcand[i]=FALSE;
}
static int intersect_row(const bool *colcand, const discrete *g1, const discrete *g2, const int ncol) {
  int i;
  int cnt = 0;
  for (i=0; i< ncol; i++)
    if (colcand[i] && (g1[i] == g2[i]) && (g1[i]!=0)) cnt++;
  return cnt;
}
static int reverse_row(const bool *colcand, const discrete *g1, const discrete *g2,
		       const int ncol, const int *symbols) {
  int i;
  int cnt = 0;
  for (i = 0; i < ncol; i++) {
    if (colcand[i] && (symbols[g1[i]] == -symbols[g2[i]])) cnt++;
  }
  return cnt;
} 

/* update seed */
void seed_update (const discrete *s, const int ncol, bits16 **profile) {
  int i;
  for(i=0;i<ncol;i++) profile[i][s[i]]++;
}

/* calculate the coverage of any row to the current consensus
 * cnt = # of valid consensus columns
 */
static void seed_current_modify (const discrete *s, bool *colcand, 
				 int* cnt, int components,
				 const int ncol, const int sigma, const double tolerance, bits16** profile) {
	int i, j, flag, n;
	int threshold = ceil(components * tolerance);
	discrete ss;
	*cnt = 0;

	for (i=0; i<ncol; i++) {
	  flag = 0; ss = s[i];
	  for (j=1; j<sigma; j++) {
	    n = profile[i][j];
	    if (j == ss) n++;
	    if (n >= threshold){
	      flag = j; break;
	    }
	  }
	  if (flag) {
	    (*cnt)++;
	    colcand[i] = TRUE;
	  }
	}
}

static bool check_seed(Edge *e, Block **bb, 
		       const int block_id,
		       const int nrow) {
  int profs[nrow];
  int i,b1,b2,b3;
  bool flag = FALSE;
  b1 = b2 = -1;
  
  /* case 1: both genes in the stack */
  for (i = 0; i < block_id; i++)
    if ( isBothInStack(bb[i]->genes, e->gene_one, e->gene_two) ) 
      return FALSE; 
  
  /* case 2: both genes NOT in the stack */
  /* iterate over blocks and find gene_one*/
  for ( i = 0; i < block_id; i++)
    if ( isInStack(bb[i]->genes, e->gene_one) ) {flag = TRUE; break;}
  if (flag) b1 = i;

  /* iterate over blocks and find gene_two*/
  flag = FALSE;	
  for ( i = 0; i < block_id; i++)
    if ( isInStack(bb[i]->genes, e->gene_two) ) {flag = TRUE; break;}
  if (flag) b2 = i;

  /* if neither one was seen */
  if ( (b1 == -1)||(b2 == -1) ) 
    return TRUE;
  /* both seen */
  /* case 3: only one gene in the stack */
  else {
      for ( i = 0; i < nrow; i++) profs[i] = 0;

      for ( i = 0; i < bb[b1]->block_rows; i++)
	profs[dsItem(bb[b1]->genes,i)]++;
      for ( i = 0; i < bb[b2]->block_rows; i++)
	profs[dsItem(bb[b2]->genes,i)]++;
      
      for ( i = 0; i < nrow; i++)
	if (profs[i] > 1) return FALSE;
      b3 = MAX(bb[b1]->block_cols, bb[b2]->block_cols);
      if ( e->score <b3/* (bb[b1]->block_cols + bb[b2]->block_cols) / 2*/ ) return FALSE;
      else return TRUE;
    }
  error("[Error] never see this message\n");
  return FALSE;
}

static void block_init(Edge *e, Block *b, 
		       dyStack *genes, dyStack *scores,
		       bool *candidates, const int cand_threshold,
		       int *components, dyStack *allincluster,
		       const int nrow, const int ncol, 
		       discrete ** const arr_c, const int mincol) {
  int i,score,top;
  int cnt = 0;
  int max_cnt, max_i;
  int *arr_rows, *arr_rows_b;
  AllocArray(arr_rows, nrow);
  AllocArray(arr_rows_b, nrow);	
  bool *colcand;
  AllocArray(colcand, ncol);
  for (i=0; i< ncol; i++) colcand[i] = FALSE;

  discrete *g1, *g2;
  g1 = arr_c[dsItem(genes,0)];
  g2 = arr_c[dsItem(genes,1)];

  /* edge score of the seed */
  for (i=0; i< ncol; i++)
    if ((g1[i]==g2[i])&&(g1[i]!=0)) colcand[i] = TRUE;
  
  /* arr_rows (arr_rows_b having the same content) records the number of sharing conditions with gene1 */
  for (i = 0; i < nrow; i++) {
      arr_rows[i] = intersect_row(colcand, g1, arr_c[i], ncol);
      arr_rows_b[i] = arr_rows[i];
  }

  /* in case more than 100 rows have the same patterns, only pick the first 100
   * , note that qsort sorts the items in an ascending order
  */
  if (nrow > 100) {
    qsort(arr_rows_b, nrow, sizeof *arr_rows, compare_int);
    top = arr_rows_b[nrow-100];
    for (i = 0; i < nrow; i++)
      if (arr_rows[i] < top) candidates[i] = FALSE;
  }
  
  while (*components < nrow) {
    max_cnt = -1;
    max_i = -1;
    (*components)++;
    for (i=0; i< nrow; i++) {
      if (!candidates[i]) continue;
	
      cnt = intersect_row(colcand,g1,arr_c[i], ncol);
      if (cnt < cand_threshold) candidates[i] = FALSE;
      if (cnt > max_cnt) {
	max_cnt = cnt;
	max_i = i;
      }
    }
    
    if (max_cnt < mincol || max_i < 0) {
      break;
    } else {	
      score = MIN(*components, max_cnt);
      if (score > b->score) b->score = score;
      dsPush(genes, max_i);
      dsPush(scores,score);
      update_colcand(colcand,g1, arr_c[max_i], ncol);
      candidates[max_i] = FALSE;
    }
  }
  free(colcand);
  free(arr_rows);
  free(arr_rows_b);
}

/* stack functions */

/* Initialize a stack */
struct dyStack *dsNew(const int size) {
  int stackSize = (size+1) * sizeof(int);
  struct dyStack *ds = malloc(stackSize);
  dsClear(ds);
  return ds;
}

/* push operation on a stack */
void dsPush(struct dyStack *ds, int element) {
    ds->items[++ds->top] = element;
}

/* Return the number of common components between two dynamic stacks */
int dsIntersect(struct dyStack *ds1, struct dyStack *ds2) {
  int cnt = 0;
  int i;
  
  for (i=0; i<dsSize(ds1); i++)
    if (isInStack(ds2, ds1->items[i])) cnt++;
  
  return cnt;
}


/* compare function for qsort, descending by score */
static int compare_block(const void *a, const void *b) {
  return ((*(Block **)b)->score - (*(Block **)a)->score);
}

static void sort_block_list(Block **el, int n) {
  qsort(el, n, sizeof *el, compare_block);
}

/* Test whether an elemente is in stack */
bool isInStack(struct dyStack *ds, int element) {
    bool flag = FALSE;
    int i;
    for (i=0; i<dsSize(ds); i++) {
      if (ds->items[i]==element) {
	  flag = TRUE; break;
	}
    }
    return flag;
}

bool isBothInStack(struct dyStack *ds, int element1, int element2) {
  bool flag1 = FALSE;
  bool flag2 = FALSE;
  bool flag = FALSE;
  int i;
  for (i=0;i<dsSize(ds); i++) {
    if (ds->items[i]==element1)
      flag1 = TRUE; 
    if (ds->items[i]==element2)
      flag2 = TRUE;
    if(flag1 && flag2) {
      flag = TRUE; break;
    }
  }
  return(flag);
}

void scan_block (struct dyStack *gene_set, Block *b_ptr,
		 discrete **arr_c,
		 const int ncol,
		 const int sigma, const double tolerance, bits16 **profile) {
  int i, j;
  
  int block_rows, cur_rows;
  block_rows = cur_rows = dsSize(gene_set);
  
  int k;
  for (j = 0; j < ncol; j++)
    for (k=0; k<sigma; k++) 
      profile[j][k] = 0;
  for (j = 0; j< cur_rows ; j++)
    seed_update(arr_c[dsItem(gene_set,j)], ncol, profile);

  int btolerance = ceil(tolerance* block_rows);
  for (j = 0; j < ncol; j++) {
    /* See if this column satisfies tolerance */
    for (i = 1; i < sigma; i++) {
      /* Pay attention to the ignored char '.' */
      if ((profile[j][i] >= btolerance)/* && (i != ignore_index) && (arr_c[dsItem(gene_set,0)][j] == i)*/) {
	dsPush(b_ptr->conds, j); break;
      }
    }		
  }
  
  b_ptr->block_cols = dsSize(b_ptr->conds);
}

/**************************************************************/
/* seed functions */
/* remove a row from the profile */
void seed_deduct (const discrete *s, const int ncol, discrete **profile) {
  int i;
  for(i=0; i<ncol; i++) profile[i][s[i]]--;
}

void seed_intersect (const discrete *s1, const discrete *s2,
		     const int ncol, bits16 **profile) {
  seed_update(s1, ncol, profile);
  seed_update(s2, ncol, profile);
}

/* coverage of any row to the current consensus 
 * cnt = number of valid consensus columns
*/
void seed_current (const discrete *s, int* cnt, int components, 
		   const int ncol, const int sigma, bits16 **profile) {
  int i, j, flag, n;
  int threshold = components;
  discrete ss;
  *cnt=0;

  for(i=0;i<ncol;i++) {
    flag=0; ss=s[i];
    for(j=1; j<sigma; j++) {
      n = profile[i][j];
      if (j==ss) n++;
      if(n>threshold) {
	flag=j; break;
      }
    }
    if (flag) (*cnt)++;
  }
}

/* coverage of any row to the current consensus
 * cnt = number of valid consus columns
 * m_cnt = number of hits of current row to the consensus
 * Note that m_cnt <= cnt always hold true
 */
void seed_intersect_r (const discrete *s, int* cnt, int* m_cnt,
		       int components, 
		       const int ncol, const int sigma, const double tolerance,
		       const bits16 **profile) {
  int i, j, flag, n;
  int threshold = ceil(components * tolerance);
  discrete ss;
  *m_cnt = *cnt = 0;
  
  for(i=0;i<ncol;i++) {
    flag=0; ss=s[i];
    for(j=0;j<sigma;j++) {
      n=profile[i][j];
      if(j==ss) n++;
      if(n>=threshold) {
	flag=j; break;
      }
    }
    if(flag) {
      (*cnt)++;
      if (flag=ss) (*m_cnt)++;
    }
  }
}



/**************************************************************/
/* core bicluster algorithm */
/* int cluster (Edge **el, int n, 
 *	     const int report_block_no,  
 *	     const discrete **dexp, const int nrow, const int ncol, 
 *	     const int sigma, 
 *	     const int mincol, const double tolerance, 
 *	     const double filter_proportion) {
 */
SEXP qubicluster (SEXP seeds, SEXP exprs, SEXP sigma_val, SEXP symbols_val, SEXP report_no,
	      SEXP tolerance_val, SEXP filter_proportion) {

  int i, j, k, ki, components;

  /*----------handle input parameters----------*/

  /* unwire the seeds */
  /* R-note: Note that here one must first extract the attribute to get the external address 
   * otherwise the returned pointer is a wrong one
   */
  ppEdge el = R_ExternalPtrAddr(getAttrib(seeds, install(EDGELIST_ATT_NAME)));
  if(el == NULL) {
    error("seeds is not a legal object\n");
    return 0;
  }
  int n = INTEGER(seeds)[0];
  int mincol = INTEGER(getAttrib(seeds, install(MINCOL_ATT_NAME)))[0];
  
  /* quantileDiscretized expression matrix */
  int *dexprs;
  discrete **dexp;
  
  /* dim and values of the matrix */
  SEXP dim;
  PROTECT(dim = getAttrib(exprs, R_DimSymbol));
  int nrow = INTEGER(dim)[0];
  int ncol = INTEGER(dim)[1];
  UNPROTECT(1);
  
  dexprs = INTEGER(exprs);
  AllocArray(dexp, nrow);
  for(i=0; i<nrow;i++) {
    AllocArray(dexp[i], ncol);
    for(j=0;j<ncol;j++) {
      dexp[i][j] = dexprs[i + nrow * j];
    }
  }
  
  /* other parameters*/
  int sigma = INTEGER(sigma_val)[0];
  int *symbols;
  int *dsymbol = INTEGER(symbols_val);
  AllocArray(symbols, sigma);
  for(i=0; i<sigma; i++)
    symbols[i] = dsymbol[i];

  int reportBlockNo = INTEGER(report_no)[0];
  double tolerance = REAL(tolerance_val)[0];
  double filterProportion = REAL(filter_proportion)[0];

  /*----------data structures----------*/
  /* initialize the block counter */
  int block_id = 0;

  /* array of block pointers
   * double size of the report block number
   */
  Block **bb;
  int allocated = reportBlockNo * 2;
  AllocArray(bb, allocated); 
  
  pEdge e;
  Block *b;
  dyStack *genes, *scores, *b_genes, *allincluster;
  bits16 **profile;

  /* profile: ncol*sigma 2-way array */
  AllocArray(profile, ncol);
  for (j = 0; j < ncol; j++) AllocArray(profile[j], sigma);
  
  genes = dsNew(nrow);
  scores = dsNew(nrow);
  allincluster = dsNew(nrow);
  
  bool *candidates;
  AllocArray(candidates, nrow);
  
  /*----------builds up blocks----------*/
  e=*el; i = 0;
  while (i++ < n) {
    e = *el++;
    /* check if both genes already enumerated in previous blocks */
    bool flag = TRUE;

    /* speed up the program if the nrow bigger than 200. WHY? */
    if (nrow > 200) {
      if ( isBothInStack(allincluster, e->gene_one, e->gene_two) )
	flag = FALSE;
    } else {
      flag = check_seed(e, bb, block_id, nrow);
    }
    if (!flag) continue;
    
    for (j = 0; j < ncol; j++)
      for (k = 0; k < sigma; k++) 
	profile[j][k] = 0;
    
    AllocVar(b);
    b->score = MIN(2, e->score);
    
    /* initialize the stacks genes and scores */
    /* set stack pointer to -1 */
    dsClear(genes);
    dsClear(scores);
    /* fill the stack */
    for(j = 0; j < nrow; j++) {
      dsPush(genes,-1);
      dsPush(scores,-1);
    }
    /* set stack pointer to -1 */
    dsClear(genes);
    dsClear(scores);
    
    dsPush(genes, e->gene_one);
    dsPush(genes, e->gene_two);
    dsPush(scores, 1);
    dsPush(scores, b->score);
    
    /* branch-and-cut condition for seed expansion */
    int cand_threshold = floor(mincol * tolerance);
    if (cand_threshold < 2) cand_threshold = 2;
    
    /* maintain a candidate list to avoid looping through all rows */		
    for (j = 0; j < nrow; j++) candidates[j] = TRUE;
    candidates[e->gene_one] = candidates[e->gene_two] = FALSE;
    
    components = 2;
    
    /* expansion step, generate a bicluster without noise */
    block_init(e, b, genes, scores, candidates, 
	       cand_threshold, &components, allincluster,
	       nrow, ncol, dexp, mincol);
    
    /* track back to find the best score that which genes makes it */
    for(k = 0; k < components; k++)
      if ((dsItem(scores,k) == b->score)&&(dsItem(scores,k+1)!= b->score)) 
	break;
    components = k + 1;
    
    for (ki=0; ki < nrow; ki++)
      candidates[ki] = TRUE;
    
    for (ki=0; ki < components - 1 ; ki++) {
      seed_update(dexp[dsItem(genes,ki)], ncol, profile);
      candidates[dsItem(genes,ki)] = FALSE;
    }
    candidates[dsItem(genes,k)] = FALSE;
    genes->top = k ;
    int cnt = 0;
    bool *colcand;
    AllocArray(colcand, ncol);
    for(ki = 0; ki < ncol; ki++) colcand[ki] = FALSE;             
    
    /* add columns satisfy the conservative r */ 
    seed_current_modify(dexp[dsItem(genes,k)], colcand, 
			&cnt, components, 
			ncol, sigma, tolerance, profile);
    
    /* add some new possible genes */
    int m_cnt;
    int tol_cnt = floor(cnt * tolerance);
    for ( ki = 0; ki < nrow; ki++) {
      m_cnt = intersect_row(colcand, dexp[dsItem(genes,0)], dexp[ki], ncol);
      if ( candidates[ki] && (m_cnt >= tol_cnt) ) {
	dsPush(genes,ki);
	components++;
	candidates[ki] = FALSE;
      }
    }
    b->block_rows_pre = components;
    
    /* add genes that negative regulated to the consensus */
    for ( ki = 0; ki < nrow; ki++) {
      m_cnt = reverse_row(colcand, dexp[dsItem(genes,0)], dexp[ki], ncol, symbols);
      if ( candidates[ki] && (m_cnt >= tol_cnt) ) {
	dsPush(genes,ki);
	components++;
	candidates[ki] = FALSE;
      }
    }
    free(colcand);
      
    /* save the current cluster*/
    b_genes = dsNew(b->block_rows_pre);
    for (ki = 0; ki < b->block_rows_pre; ki++)
      dsPush(b_genes, dsItem(genes,ki));
    
    /* store gene arrays inside block */
    b->genes = dsNew(components);
    b->conds = dsNew(ncol);
    
    /* add tolerated columns */
    scan_block(b_genes, b, 
	       dexp, ncol, 
	       sigma, tolerance, profile);
    if (b->block_cols == 0) continue;

    b->block_rows = components;
    b->score = b->block_rows * b->block_cols;		
    dsClear(b->genes);
    for ( ki=0; ki < components; ki++)
      dsPush(b->genes,dsItem(genes,ki));
    for(ki = 0; ki < components; ki++)
      if(!isInStack(allincluster, dsItem(genes,ki))) 
	dsPush(allincluster,dsItem(genes,ki));	
    
    bb[block_id++] = b;
    
    /* reaching the results number limit */
    if (block_id == allocated) break;
  }
  
  // TODO: return the index of rows/columns of each cluster
  SEXP rb=report_blocks(bb, block_id, reportBlockNo, filterProportion, dexp, symbols);

  /* free-up the candidate list */
  free(bb);
  free(candidates);
  free(allincluster);
  free(dexp);
  free(profile);
  free(symbols);
  
  return(rb);
}

SEXP report_blocks(Block** bb, 
		  const int num, 
		  const int reportBlockNo,
		  const double filter_proportion,
		  discrete **arr_c, const int *symbols) {
  /* sort blocks by scores */
  sort_block_list(bb, num);
  
  int i,j,k, bc, p;
  int n = MIN(num, reportBlockNo);
  bool flag;
  
  Block **output;
  AllocArray(output, n);
  
  Block **bb_ptr = output;
  Block *b_ptr;
  double cur_rows, cur_cols;
  double inter_rows, inter_cols;

  int **blocks_rows;
  int *blocks_rows_cnt;
  int **blocks_cols;
  int *blocks_cols_cnt;
  AllocArray(blocks_rows_cnt, n);
  AllocArray(blocks_cols_cnt, n);
  AllocArray(blocks_rows, n);
  AllocArray(blocks_cols, n);

  /* Filter overlapping blocks */
  i = j = 0;
  while (i < num && j < n) {
    b_ptr = bb[i];
    cur_rows = b_ptr->block_rows;
    cur_cols = b_ptr->block_cols;
    
    flag = TRUE;
    
    k = 0;
    while (k < j) {
      inter_rows = dsIntersect(output[k]->genes, b_ptr->genes);
      inter_cols = dsIntersect(output[k]->conds, b_ptr->conds);
      
      if (inter_rows*inter_cols > filter_proportion*cur_rows*cur_cols) {
	flag = FALSE; break;
      }
      /*proportion=(inter_rows*inter_cols)/(cur_rows*cur_cols);
	printf ("%d\t%d\t%.3f\n",j,k,proportion);*/
      k++;
    }
    i++;

    if (flag) {
      blocks_rows_cnt[j] = b_ptr->block_rows;
      blocks_cols_cnt[j] = b_ptr->block_cols;
      AllocArray(blocks_rows[j], blocks_rows_cnt[j]);
      AllocArray(blocks_cols[j], blocks_cols_cnt[j]);
      for(bc=0; bc<dsSize(b_ptr->genes); bc++)
	blocks_rows[j][bc]=dsItem(b_ptr->genes, bc);
      for(bc=0; bc<dsSize(b_ptr->conds); bc++)
	blocks_cols[j][bc]=dsItem(b_ptr->conds, bc);
      j++;
      *bb_ptr++ = b_ptr;
    }
  }

  SEXP ans, rowIndex, colIndex, ansItem;
  
  PROTECT(ans=allocVector(VECSXP, j));
  for(bc=0; bc<j; bc++) {
    PROTECT(rowIndex=allocVector(INTSXP, blocks_rows_cnt[bc]));
    PROTECT(colIndex=allocVector(INTSXP, blocks_cols_cnt[bc]));
    /* plus one since the R index begins from 1*/
    for(p=0;p<blocks_rows_cnt[bc];p++)
      INTEGER(rowIndex)[p]=blocks_rows[bc][p]+1;
    for(p=0;p<blocks_cols_cnt[bc];p++)
      INTEGER(colIndex)[p]=blocks_cols[bc][p]+1;
    PROTECT(ansItem=allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ansItem, 0, rowIndex);
    SET_VECTOR_ELT(ansItem, 1, colIndex);
    SET_VECTOR_ELT(ans, bc, ansItem);
    UNPROTECT(3);
  }
  UNPROTECT(1);
  return(ans);
}
