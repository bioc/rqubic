/* quantile_from_sorted_data, modified from the function quantile_from_sorted_data
 * implemented in the GNU Scientific Library.
*/

#include "rqubic.h"
#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}
#define REGISTER_CCALLABLE(fun) \
        R_RegisterCCallable("rqubic", #fun, (DL_FUNC) &fun)


/* Safe wrapper for standard malloc */
void* xmalloc ( int size )
{
  register void* value = malloc(size);
  if (value == NULL)
    REprintf("[Error] Memory exhausted (xmalloc)");
  return value;
}

/* Safe wrapper for standard realloc */
void* xrealloc ( void* ptr, int size )
{
  register void* value = realloc(ptr, size);
  if (value == NULL)
    REprintf("[Error] Memory exhausted (xrealloc)");
  return value;
}

static const R_CallMethodDef callMethods[] = {

  /* discrete.c */
  CALLMETHOD_DEF(discretize_matrix,3),

  /* graph.c */
  CALLMETHOD_DEF(generate_sorted_seeds,2),

  /* bicluster.c */
  CALLMETHOD_DEF(qubicluster, 7),

  {NULL, NULL, 0}
};

static const R_CMethodDef cMethods[] ={
  {NULL, NULL, 0}
};

void R_init_rqubic(DllInfo *info)
{
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
}


/* initializer */
SEXP RQUBIC_init(SEXP env) {
  RQUBIC_edgelist_tag = install("RQUBIC_edgelist");
  return R_NilValue;
}
