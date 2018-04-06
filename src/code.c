/* History: Jan 26 2011 Make the breaking of 3+ ties more efficient
            Feb 23 2012 Add code for the new method of generating permutations
            Mar 09 2012 Add code to only store the snps in each gene defined by inspect.vec.gene.list
                        This will reduce the amount of memory needed.
            Jul 30 2013 Remove printf/exit calls for CRAN submission
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <Rinternals.h>

#define CHECK_MEM(obj) if (obj == NULL) {error("ERROR: allocating memory \n");}
#define MAX( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define MIN( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define SMALL_DOUBLE -1.0e300
#define LARGE_DOUBLE 1.0e300
#define SMALL_POSITIVE 1.0e-300
#define TIES_EPS 1.0e-15
#define SMALL_INT -1e8
#define DEBUG 0

typedef unsigned char binary;

void ARTP_pathway(char **obsfile, char **permfile, int *nr, int *nc, int *cols, int *ncols, int *geneStart, int *geneStop, int *ngene, int *rowsize,\
                  int *inspect_list, int *llen, int *iStart, int *iStop, int *inspect_path, int *ilen, char **r2Files, int *r2Flag, char **outfile,\
                  int *ties_method, int *method, double *parms, double *ret_p_obs, double *ret_minp, int *ret_nperm);



static const R_CMethodDef callMethods[] = {
  {"ARTP_pathway", (DL_FUNC)&ARTP_pathway, 25},
  {NULL, NULL, 0}
};




/*
double dvec_max(double *, int, int *);
void ivecinit(int *, int, int);
int ivecsum(int *, int);
void bvecinit(binary *, int, binary);
int bvecsum(binary *, int);
*/

/**********************************************************/
/* Functions for vectors */
/**********************************************************/

/* Function to sum a binary vector */
static int bvecsum(vec, n)
binary *vec;
int n;
{
  int sum = 0;
  
  while (n-- > 0) {
    sum += *vec++;
  } 
  return(sum);

}

/* Function to check the bounds of a vector. Changes elements if too small */
static int checkVecBounds(vec, n, lower, upper)
double *vec;
int n;
double lower;
double upper;
{
  int i;
  double *ptrd, temp;
  
  for (i=0, ptrd=vec; i<n; i++, ptrd++) {
    temp = *ptrd;
    if ((temp < lower) || (temp > upper)) return(1); 
    if (temp < SMALL_POSITIVE) *ptrd = SMALL_POSITIVE;
  }
  return(0);
}

/* Function to initialize a binary vector to a constant */
static void bvecinit(vec, n, value)
binary *vec, value;
int n;
{
  binary *ptr = vec;
  while (n-- > 0) {
    *ptr++ = value;
  }
}

/* Function to initialize a double vector to another vector */
static void dvec_dvec_init(outvec, vec, n)
double *outvec, *vec;
int n;
{
  double *ptr1 = vec, *ptr2 = outvec;
  while (n-- > 0) {
    *ptr2++ = *ptr1++;
  }
}

/* Function to reverse a double vector */
static void dvec_reverse(outvec, vec, n)
double *outvec, *vec;
int n;
{
  int i, nm1;
  double *ptr, *p2;
  
  nm1 = n-1;
  /*for (i=0, ptr=outvec; i<n; i++, ptr++) *ptr = vec[nm1-i];*/
  for (i=0, ptr=outvec, p2=vec+nm1; i<n; i++, ptr++, p2--) *ptr = *p2;

}

/* Function to get the minimum of a double vector */
static double dvec_min(vec, n)
double *vec;
int n;
{
  double *ptr, minval=LARGE_DOUBLE, temp;
  int i;

  for (i=0, ptr=vec; i<n; i++, ptr++) {
    temp = *ptr;
    if (temp < minval) minval = temp;
  }

  return(minval);
}

/* Function to get the maximum and index of the maximum of a double vector */
static double dvec_max(vec, n, index)
double *vec;
int n, *index;
{
  double *ptr, maxval=SMALL_DOUBLE, temp;
  int i, jj=0;

  for (i=0, ptr=vec; i<n; i++, ptr++) {
    temp = *ptr;
    if (temp > maxval) {
      maxval = temp;
      jj     = i;
    }
  }
  *index = jj;
  return(maxval);
}

/* Function to get the maximum and index of the maximum of an integer vector */
static int ivec_max(vec, n, index)
int *vec;
int n, *index;
{
  int *ptr, maxval=SMALL_INT, temp;
  int i, jj=0;

  for (i=0, ptr=vec; i<n; i++, ptr++) {
    temp = *ptr;
    if (temp > maxval) {
      maxval = temp;
      jj     = i;
    }
  }
  *index = jj;
  return(maxval);
}

/* Function to compute the cumulative sum of a vector */
static void cumsum(x, n, ret)
double *x, *ret;
int n;
{
  int i;
  double sum=0.0, *ptrd, *ptr;

  for (i=0, ptrd=x, ptr=ret; i<n; i++, ptrd++, ptr++) {
    sum += *ptrd;
    *ptr = sum;
  }

} /* END: cumsum */


/*
void print_iVec(vec, n, name)
int *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %d ", vec[i]);
  }
  printf("\n \n");
}

void print_dVec(vec, n, name)
double *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %g ", vec[i]);
  }
  printf("\n \n");
}

void print_dMat(mat, nr, nc, name)
double **mat;
int nr, nc;
char name[10];
{
  int i, j;
  printf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) printf(" %g ", mat[i][j]);
    printf("\n");
  }
  printf("\n \n");
}

void print_iMat(mat, nr, nc, name)
int **mat;
int nr, nc;
char name[10];
{
  int i, j;
  printf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) printf(" %d ", mat[i][j]);
    printf("\n");
  }
  printf("\n \n");
}

void print_bMat(mat, nr, nc, name)
binary **mat;
int nr, nc;
char name[10];
{
  int i, j;
  printf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) printf(" %d ", mat[i][j]);
    printf("\n");
  }
  printf("\n \n");
}
*/

/* Sorting function, x is changed, sort from small to large*/
static void quicksort(int start, int stop, double *x, int *cvec)
{
 int i, j, k;
 double temp, median;
 int tempd;


  while (start < stop) {
    /*
    ** first-- if the list is short, do an ordinary insertion sort
    */
    if ((stop-start)<11) {
	for (i=start+1; i<=stop; i++) {
	    temp = x[i];
	    tempd= cvec[i];
	    j=i-1;

	    while (j>=start && (x[j]>temp)) {
		x[j+1] = x[j];
		cvec[j+1] = cvec[j];
		j--;
		}
	    x[j+1] = temp;
	    cvec[j+1]  = tempd;
	    }
	return;
	}

    /*
    ** list is longer -- split it into two
    **  I use the median of 3 values as the split point
    */
    i=start;
    j=stop;
    k = (start + stop)/2;

    median = x[k];
    if (x[i] >= x[k]) {      /* one of j or k is smallest */
	if (x[j] > x[k]) {   /* k is smallest */
	    if (x[i] > x[j])  median = x[j];
	    else median= x[i];
	    }
	}
    else {
	if (x[j] < x[k]) {
	    if (x[i] > x[j]) median = x[i];
	    else median = x[j];
	    }
	}

    /* 
    **  Now actually do the partitioning 
    **   Because we must have at least one element >= median, "i"
    **   will never run over the end of the array.  Similar logic
    **   applies to j.
    ** A note on the use of "<" rather than "<=".  If a list has lots
    **   of identical elements, e.g. 80/100 are "3.5", then we will
    **   often go to the swap step with x[i]=x[j]=median.  But we will
    **   get the pointers i and j to meet approximately in the middle of
    **   the list, and that is THE important condition for speed in a
    **   quicksort.
    **   
    */
    while (i<j) {
	/*
	** top pointer down till it points at something too large
	*/
	while (x[i] < median) i++;

	/*
	** bottom pointer up until it points at something too small
	*/
	while(x[j] > median) j--;

	if (i<j) {
	    if (x[i] > x[j]) {  /* swap */
		temp = x[i];
		x[i] = x[j];
		x[j] = temp;
		tempd= cvec[i];   cvec[i] =cvec[j];  cvec[j] =tempd;
		}
	    i++; j--;
	    }
	}

    /*
    ** The while() step helps if there are lots of ties.  It will break
    **  the list into 3 parts: < median, ==median, >=median, of which only
    **  the top and bottom ones need further attention.
    ** The ">=" is needed because i may be  == to j
    */
    while (x[i] >= median && i>start) i--;
    while (x[j] <= median && j<stop ) j++;

    /*
    ** list has been split, now do a recursive call
    **   always recur on the shorter list, as this keeps the total
    **       depth of nested calls to less than log_base2(n).
    */
    if ((i-start) < (stop-j)) { /* top list is shorter */
	if ((i-start)>0) quicksort(start,i, x, cvec);
	start =j; 
	}

    else {    /* bottom list is shorter */
	if ((stop -j)>0) quicksort(j,stop, x, cvec);
	stop=i; 
	}
     }
}

/* Sorting function for an in teger vector, x is changed, sort from small to large*/
static void iquicksort(int start, int stop, int *x, int *cvec)
{
 int i, j, k;
 int temp, median;
 int tempd;


  while (start < stop) {
    /*
    ** first-- if the list is short, do an ordinary insertion sort
    */
    if ((stop-start)<11) {
	for (i=start+1; i<=stop; i++) {
	    temp = x[i];
	    tempd= cvec[i];
	    j=i-1;

	    while (j>=start && (x[j]>temp)) {
		x[j+1] = x[j];
		cvec[j+1] = cvec[j];
		j--;
		}
	    x[j+1] = temp;
	    cvec[j+1]  = tempd;
	    }
	return;
	}

    /*
    ** list is longer -- split it into two
    **  I use the median of 3 values as the split point
    */
    i=start;
    j=stop;
    k = (start + stop)/2;

    median = x[k];
    if (x[i] >= x[k]) {      /* one of j or k is smallest */
	if (x[j] > x[k]) {   /* k is smallest */
	    if (x[i] > x[j])  median = x[j];
	    else median= x[i];
	    }
	}
    else {
	if (x[j] < x[k]) {
	    if (x[i] > x[j]) median = x[i];
	    else median = x[j];
	    }
	}

    /* 
    **  Now actually do the partitioning 
    **   Because we must have at least one element >= median, "i"
    **   will never run over the end of the array.  Similar logic
    **   applies to j.
    ** A note on the use of "<" rather than "<=".  If a list has lots
    **   of identical elements, e.g. 80/100 are "3.5", then we will
    **   often go to the swap step with x[i]=x[j]=median.  But we will
    **   get the pointers i and j to meet approximately in the middle of
    **   the list, and that is THE important condition for speed in a
    **   quicksort.
    **   
    */
    while (i<j) {
	/*
	** top pointer down till it points at something too large
	*/
	while (x[i] < median) i++;

	/*
	** bottom pointer up until it points at something too small
	*/
	while(x[j] > median) j--;

	if (i<j) {
	    if (x[i] > x[j]) {  /* swap */
		temp = x[i];
		x[i] = x[j];
		x[j] = temp;
		tempd= cvec[i];   cvec[i] =cvec[j];  cvec[j] =tempd;
		}
	    i++; j--;
	    }
	}

    /*
    ** The while() step helps if there are lots of ties.  It will break
    **  the list into 3 parts: < median, ==median, >=median, of which only
    **  the top and bottom ones need further attention.
    ** The ">=" is needed because i may be  == to j
    */
    while (x[i] >= median && i>start) i--;
    while (x[j] <= median && j<stop ) j++;

    /*
    ** list has been split, now do a recursive call
    **   always recur on the shorter list, as this keeps the total
    **       depth of nested calls to less than log_base2(n).
    */
    if ((i-start) < (stop-j)) { /* top list is shorter */
	if ((i-start)>0) iquicksort(start,i, x, cvec);
	start =j; 
	}

    else {    /* bottom list is shorter */
	if ((stop -j)>0) iquicksort(j,stop, x, cvec);
	stop=i; 
	}
     }

} /* END: iquicksort */ 

/*******************************************************/
/* Functions to break ties from a vector of ranks      */
/*******************************************************/
static void breakTies(rnk, n, method)
int *rnk, n, method;
{
  int i, *tempi, rnum, temp; 

  if (!method) {
    /* random */
    if (n == 2) {
      /* Generate a uniform random deviate and switch if < 0.5 */
      if (unif_rand() < 0.5) {
        i = *rnk;
        *rnk = rnk[1];
        rnk[1] = i;
      }
    } else {
      /* Generate a uniform number between 0 and n-1, and swap with ith position */
      for (i=0; i<n; i++) {
        rnum = floor(n*unif_rand());
        if (rnum == n) rnum = rnum - 1;
        temp = rnk[i];
        rnk[i] = rnk[rnum];
        rnk[rnum] = temp;
      }
    }
  } else {
    tempi = (int *) malloc(n*sizeof(int));
    CHECK_MEM(tempi);

    /* First */
    iquicksort(0, n-1, rnk, tempi);

    free(tempi);
  }

} /* END: breakTies */

static void ties(vec, ivec, n, method)
double *vec; /* Must be sorted */
int *ivec, n, method;
{
  int i, begin=0, end=0, flag, len, *pint=0;
  double *p0, *p1;

  if (n == 1) return;

  flag = 0;
  p0   = vec;
  p1   = vec + 1;
  for (i=0; i<n-1; i++) {
    if (fabs(*p0 - *p1) < TIES_EPS) {
      if (!flag) {
        flag  = 1;
        begin = i; 
        pint  = &ivec[begin];
      } 
      end = i + 1;
    } else if (flag) {
      /* Break the ties from begin to end */
      len = end - begin + 1;
      breakTies(pint, len, method);
      flag = 0;
    }
    p0++;
    p1++;
  }
  if (flag) {
    /* Break the ties from begin to end */
    len = end - begin + 1;
    breakTies(pint, len, method);
  }

} /* END: breakTies */

/* Function to get the ranks of a double vector */
static void get_ranks(vec, n, rank, ties_method)
double *vec;
int n;
int *rank; /* Output ranks */
int ties_method; /* 1 = first, 2 = random */
{
  int i, *pi1, *tempi;
  double *temp, *p1, *p2;

  tempi = (int *) malloc(n*sizeof(int));
  CHECK_MEM(tempi);
  temp = (double *) malloc(n*sizeof(double));
  CHECK_MEM(temp);

  /* Copy vec to temp */
  for (i=0, p1=vec, p2=temp, pi1=tempi; i<n; i++, p1++, p2++, pi1++) {
    *p2 = *p1;
    *pi1 = i;
  }

  quicksort(0, n-1, temp, tempi);

  /* Break the ties */
  ties(temp, tempi, n, ties_method);
  free(temp);

  /* Get the correct order for the vector */
  for (i=0, pi1=tempi; i<n; i++, pi1++) rank[*pi1] = i;
  free(tempi);
  
} 

/***********/
/* my_sort */
/***********/
static void my_sort(vec, n, r2Mat, ret_vec, accept, used, rvec, remain)
double *vec, *ret_vec;
int n;
binary **r2Mat, *accept, *used;
int *rvec;
double *remain;
{
  int i, j, k, nused, u, nremain;
  double temp, *ptrd, *p2, *ptr_ret, maxval;
  binary *ptrb, *ptrb2;

  if (n == 1) {
    *ret_vec = *vec;
    return;
  }

  /* Initialize accept to 1, used to 0 */
  bvecinit(accept, n, 1);
  bvecinit(used, n, 0);

  /* Get the snp with largest value in vec */  
  temp = dvec_max(vec, n, &u);
  nused      = 1;
  ret_vec[0] = temp;
  accept[u]  = 0;
  used[u]    = 1;
  ptr_ret    = &ret_vec[1];

  for (i=1; i<n; i++) {
    /* Get the snps with r^2 < threshold. Update accept for the current used value. */
    accept[u] = 0;
    for (k=0, ptrb2=r2Mat[u], ptrb=accept; k<n; k++, ptrb2++, ptrb++) {
      if (!*ptrb2) *ptrb = 0;
    }

    /* Get the number of snps to choose from */
    j = bvecsum(accept, n);
    
    if (!j) {
      /* No more snps to choose from that are not in high ld with the set
         of snps already chosen. So add the remaining snps by highest vec value. */
      nremain = n - nused;
      ptrd = remain;
      for (k=0, ptrb=used, p2=vec; k<n; k++, ptrb++, p2++) {
        if (!*ptrb) *ptrd++ = *p2; 
      }

      /* We need largest to smallest values */
      quicksort(0, nremain-1, remain, rvec);
      ptrd = remain + (nremain - 1); 
      for (k=0, p2=ptrd; k<nremain; k++, p2--) *ptr_ret++ = *p2;
      return;
    } 

    /* There is at least 1 snp < threshold, get the one with largest vec */
    maxval = SMALL_DOUBLE;
    u = -9999;
    for (k=0, ptrb=accept, ptrd=vec; k<n; k++, ptrb++, ptrd++) {
      if ((*ptrb) && (*ptrd > maxval)) {
        u = k;
        maxval = *ptrd;
      } 
    }
    used[u] = 1;
    nused++;
    *ptr_ret++ = maxval;
  }

} /* END: my_sort */

/*************************************************************/
/* Functions for allocating matrices */
/*************************************************************/

static binary ** bMatrix_alloc(nrow, ncol)
int nrow, ncol;
{
  binary **mat, **ptr;
  int i;
  size_t size;

  mat = (binary **) malloc(nrow*sizeof(binary *));
  CHECK_MEM(mat);
  size = ncol*sizeof(binary);
  for (i=0, ptr=mat; i<nrow; i++, ptr++) {
    *ptr = (binary *) malloc(size);
    CHECK_MEM(*ptr);
  }
  return(mat);
}

static double ** dblMatrix_alloc(nrow, ncol)
int nrow, ncol;
{
  double **mat, **ptr;
  int i;
  size_t size;

  mat = (double **) malloc(nrow*sizeof(double *));
  CHECK_MEM(mat);
  size = ncol*sizeof(double);
  for (i=0, ptr=mat; i<nrow; i++, ptr++) {
    *ptr = (double *) malloc(size);
    CHECK_MEM(*ptr);
  }
  return(mat);
}

static void bMatrix_free(mat, nrow)
binary **mat;
int nrow;
{
  int i;
  binary **ptr;
  
  for (i=0, ptr=mat; i<nrow; i++, ptr++) free(*ptr);
  free(mat);
}

static void dblMatrix_free(mat, nrow)
double **mat;
int nrow;
{
  int i;
  double **ptr;
  
  for (i=0, ptr=mat; i<nrow; i++, ptr++) free(*ptr);
  free(mat);
}

/******************************************************/
/* Functions for reading the data */
/******************************************************/

/* Function for opening a data set in read mode */
static FILE * open_file(datafile, mode)
const char *datafile, *mode;
{
  FILE *fpt;

  fpt = fopen(datafile, mode);
  if (fpt == 0) {
    error("\n Error opening file: %s \n", datafile);
    /*exit(0);*/
  }

  return(fpt);
}

/* Function for reading an r2 matrix of 0-1 */
static binary ** read_matrix(datafile, n)
const char *datafile;
int n;
{
  FILE *fpt;
  int i, j, temp, val;
  binary **mat;
  
  mat = bMatrix_alloc(n, n);

  fpt = open_file(datafile, "r");

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      temp = fscanf(fpt, "%d", &val);
      if (temp == EOF) {
        error(" \n Error with matrix file \n");
        /*exit(0);*/
      }
      mat[i][j] = (binary) val;
    }
  }
  fclose(fpt);

  return(mat);
}

/************************************************************/
/* fun_ck0: Call this function once before fun.ck is called */
/************************************************************/
static void fun_ck0(xlen, inspect_vec, ilen)
int xlen, *inspect_vec, ilen;
{
  int i, max=-9999, *ptri;

  for (i=0, ptri=inspect_vec; i<ilen; i++, ptri++) {
    if (*ptri > max) max = *ptri;
  }
  
  if (max > xlen) {
    error("\n ERROR: too large a truncation cut point \n");
    /*exit(1);*/
  }

} /* END: fun_ck0 */

/**********/
/* fun_ck */
/**********/
static void fun_ck(x, xlen, inspect_vec, ilen, ret, tempy)
double *x;
int xlen, *inspect_vec, ilen;
double *ret;  /* Vector of length ilen */
double *tempy; /* Temporary vector for efficiency */
{
  int i, *ptri;
  double *ptrd;

  /* Compute the cumulative sum of x */
  cumsum(x, xlen, tempy);

  /* Get the elements for inspect.vec */
  for (i=0, ptri=inspect_vec, ptrd=ret; i<ilen; i++, ptri++, ptrd++) *ptrd = tempy[*ptri]; 

} /* END: fun_ck */

/************/
/* snpBased */
/************/
static void snpBased(T, nr, nc, start_col, inspect_vec, ilen, r2Mat, r2Flag)
double **T;
int nr, nc, *inspect_vec, ilen;
int r2Flag, start_col;
binary **r2Mat;
{
  int i, *tempi, max_inspect;
  double *temp, **prow, *p1, *tmpd=0, *tempy;
  binary *used=0, *accept=0;

  if (DEBUG) Rprintf("Begin: snpBased\n");
  if (nc > 1) {
  
    tempi = (int *) malloc(nc*sizeof(int));
    CHECK_MEM(tempi);
    temp = (double *) malloc(nc*sizeof(double));
    CHECK_MEM(temp);

    if (r2Flag) {
      tmpd = (double *) malloc(nc*sizeof(double));
      CHECK_MEM(tmpd);
      used = (binary *) malloc(nc*sizeof(binary));
      CHECK_MEM(used);
      accept = (binary *) malloc(nc*sizeof(binary));
      CHECK_MEM(accept);
    }

    /* Sort the rows of T in decreasing order */
    for (i=0, prow=T; i<nr; i++, prow++) {
      p1 = *prow + start_col;
      dvec_dvec_init(temp, p1, nc);

      if (r2Flag) {
        my_sort(temp, nc, r2Mat, p1, accept, used, tempi, tmpd);
      } else {
        quicksort(0, nc-1, temp, tempi);
        dvec_reverse(p1, temp, nc);
      }
    }

    free(tempi);
    free(temp);
    if (r2Flag) {
      free(tmpd);
      free(used);
      free(accept);
    }
  }  
  
  /* Get the maximum value of inspect_vec, since the cumsum does not need to be for the whole vector */
  max_inspect = ivec_max(inspect_vec, ilen, &i) + 1;
  nc = MIN(nc, max_inspect);

  tempy = (double *) malloc(nc*sizeof(double));
  CHECK_MEM(tempy);

  /* Apply inspect.fun to each row of T. In fun_ck ret can be T */
  for (i=0, prow=T; i<nr; i++, prow++) {
    p1 = *prow + start_col;
    fun_ck(p1, nc, inspect_vec, ilen, p1, tempy);
  }
  free(tempy);

  if (DEBUG) Rprintf("End: snpBased\n");

} /* END: snpBased */

/**************/
/* get_pvalue */
/**************/
static void get_pvalue(dev_mat, nr, nc, start_col, ret_p_obs, ret_p_perm, ties_method)
double **dev_mat;  /* First row is for the observed */
double *ret_p_obs; /* Length nc */
double *ret_p_perm; /* Length nr */
int nr, nc, ties_method, start_col;
{
  int i, j, *rank, *ptri, ii;
  double *temp, *ptrd, one_over_nr, **p2;

  if (DEBUG) Rprintf("Begin: get_pvalue\n");

  rank = (int *) malloc(nr*sizeof(int));
  CHECK_MEM(rank);
  temp = (double *) malloc(nr*sizeof(double));
  CHECK_MEM(temp);

  one_over_nr = 1.0/nr;

  for (i=0; i<nc; i++) {
    ii = i + start_col;

    if (DEBUG) Rprintf("i=%d, nc=%d, ii=%d, start_col=%d, nr=%d, ties_method=%d\n", i, nc, ii, start_col, nr, ties_method);

    /* Get the ith column and multiply by -1 */
    for (j=0, ptrd=temp, p2=dev_mat; j<nr; j++, ptrd++, p2++) *ptrd = -*(*p2 + ii);

    /* Get the ranks, and then remember to add 1 below since C vectors start at 0, not 1 */
    get_ranks(temp, nr, rank, ties_method);

    /* Replace ith column of dev_mat with (ranks+1)/nr */
    for (j=0, ptri=rank, p2=dev_mat; j<nr; j++, ptri++, p2++) *(*p2 + ii) = (*ptri + 1)*one_over_nr;

  }

  /* Set p.obs to row 0 of dev_mat */
  
  dvec_dvec_init(ret_p_obs, &dev_mat[0][start_col], nc);

  if (nc > 1) {
    /* Apply min across each row */
    for (j=0, ptrd=temp, p2=dev_mat; j<nr; j++, ptrd++, p2++) *ptrd = dvec_min(*p2 + start_col, nc);    
    get_ranks(temp, nr, rank, ties_method);

    for (j=0, ptri=rank, ptrd=ret_p_perm; j<nr; j++, ptri++, ptrd++) *ptrd = (*ptri + 1)*one_over_nr;
  } else {
    for (j=0, ptrd=ret_p_perm, p2=dev_mat; j<nr; j++, ptrd++, p2++) *ptrd = *(*p2 + start_col);
  }

  free(temp);
  free(rank);

  if (DEBUG) Rprintf("End: get_pvalue\n");

} /* END: get_pvalue */

/***************/
/* pathway_snp */
/***************/
static void pathway_snp(T, nr, nc, start_col, inspect_vec, ilen, ret_p_obs, ret_minp, r2Mat, r2Flag, ties_method)
double **T; /* first row is for the observed data */
int nr, nc, *inspect_vec, ilen;
double *ret_p_obs, *ret_minp;
int r2Flag, ties_method, start_col; 
binary **r2Mat;
{
  double p_obs, *temp;

  if (nc > 1) {
    temp = (double *) malloc(nr*sizeof(double));
    CHECK_MEM(temp);

    snpBased(T, nr, nc, start_col, inspect_vec, ilen, r2Mat, r2Flag);
    get_pvalue(T, nr, ilen, start_col, ret_p_obs, temp, ties_method);
    *ret_minp = temp[0];
    free(temp);
  } else {
    p_obs  = exp(-T[0][start_col]);
    *ret_p_obs = p_obs;
    *ret_minp  = p_obs;
  }

} /* END: pathway_snp */

/***************/
/* my_sum_gene */
/***************/
static void my_sum_gene(T, nr, nc, start_col, inspect_vec, ilen, ret_p_obs, ret_p_perm, r2Mat, r2Flag, ties_method)
double **T; /* First row is for the observed data */
int nr, nc, *inspect_vec, ilen;
double *ret_p_obs, *ret_p_perm;
int r2Flag, ties_method;
binary **r2Mat;
int start_col; /* Column to start in matrix T */
{
  
  if (nc > 1) snpBased(T, nr, nc, start_col, inspect_vec, ilen, r2Mat, r2Flag);
  
  get_pvalue(T, nr, ilen, start_col, ret_p_obs, ret_p_perm, ties_method);

} /* END: my_sum_gene */

/*****************/
/* pathway_gene1 */
/*****************/
static void pathway_gene(T, nr, nc, inspect_list, llen, iStart, iStop, geneStart, geneStop, gs_len,\
                   inspect_path, ilen, r2Files, r2Flag, outfile, ret_p_obs, ret_minp, ties_method, method, parms) 
double **T;
int nr, nc, *inspect_list, llen, *iStart, *iStop, *geneStart, *geneStop, gs_len;
char *outfile; /* Output file for gene results */
int *inspect_path, ilen;
char **r2Files; /* Character vector of length gs_len for the r2 matrices */
double *ret_p_obs, *ret_minp;  
int r2Flag; /* 0 or 1 for the r^2 option */
int ties_method;
int method; 
double *parms; 
{
  int i, ncols, a, b, j, k, outflag, *ivec, *ptri, start_col;
  double *p_perm, *p_obs, *ptrd, **p2;  
  FILE *fpt=0;
  binary **r2Mat=0;

  i = (int) strlen(outfile);
  if ((outfile == NULL) || (!i)) {
    outflag = 0;
  } else {
    outflag = 1;
    fpt = open_file(outfile, "w");
  }

  p_perm = (double *) malloc(nr*sizeof(double));
  CHECK_MEM(p_perm);

  for (i=0; i<gs_len; i++) {
    a = geneStart[i];
    b = geneStop[i];
    ncols = b - a + 1;
    start_col = a;

    if (DEBUG) Rprintf("Gene %d, ncols = %d, start_col = %d \n", i, ncols, start_col);

    /* Get the inspect.vec for this gene */
    a = iStart[i];
    b = iStop[i];
    j = b - a + 1;
    ivec = (int *) malloc(j*sizeof(int));
    CHECK_MEM(ivec);
    for (k=a, ptri=ivec; k<=b; k++, ptri++) *ptri = inspect_list[k];
   
    p_obs = (double *) malloc(j*sizeof(double));
    CHECK_MEM(p_obs);

    /* Read the r^2 matrix */
    if (r2Flag) r2Mat = read_matrix(r2Files[i], ncols);

    if (DEBUG) Rprintf("Begin: my_sum_gene \n");
    my_sum_gene(T, nr, ncols, start_col, ivec, j, p_obs, p_perm, r2Mat, r2Flag, ties_method);
    if (DEBUG) Rprintf("End: my_sum_gene \n");

    free(ivec);
    if (r2Flag) bMatrix_free(r2Mat, ncols);

    /* Store the columns of t.gene.mat in T */
    for (k=0, ptrd=p_perm, p2=T; k<nr; k++, ptrd++, p2++) *(*p2 + i) = -log(*ptrd); 

    /* Write out gene results */
    if (outflag) {
      fprintf(fpt, "%g\n", p_perm[0]);
      for (k=0, ptrd=p_obs; k<j-1; k++, ptrd++) fprintf(fpt, "%g,", *ptrd);
      fprintf(fpt, "%g\n", *ptrd);
      for (k=0, p2=T; k<nr-1; k++, p2++) fprintf(fpt, "%g,", *(*p2 + i));
      fprintf(fpt, "%g\n", *(*p2 + i));
      k = fflush(fpt);
    }
    free(p_obs);

  } /* END: for (i=0; i<gs_len; i++) */

  if (outflag) fclose(fpt);
  free(p_perm);
  r2Flag = 0;

  /* Call pathway.snp */
  if (DEBUG) Rprintf("Begin: pathway_snp \n");
  pathway_snp(T, nr, gs_len, 0, inspect_path, ilen, ret_p_obs, ret_minp, r2Mat, r2Flag, ties_method);
  if (DEBUG) Rprintf("End: pathway_snp \n");  

} /* pathway_gene */

/* Function to read row 1 of column names */
static void read_header(fpt, nchars)
FILE *fpt;
int nchars;
{
  char line[nchars];

  fgets(line, sizeof(line), fpt);
}

/* Function to read a comma delimited line of doubles from a file */
static void parse_dline_comma(line, nc, ret_vec)
char *line;
int nc;
double *ret_vec;
{
  int i;
  double *ptrd, val;
  char *token, *end;

  ptrd = ret_vec;
  token = strtok(line, ",");
  i = 0;
 
  while (token != NULL) {
   /* strtod will return 0 if token could not be converted to a numeric value */
    val = strtod(token, &end);
    *ptrd = val;

    /* Check if value should be missing */
    if ((val < SMALL_POSITIVE) && (*end != 0) && (end == token)) *ptrd = -9999.0;
      
    ptrd++;
    i++;

    token = strtok(NULL, ",");
  }
  if (i != nc) {
    error("\n ERROR: reading data \n");  
    /*exit(1);*/
  }
}

/* Function to read the data */
static int read_data(obsfile, permfile, nr, nc, cols, ncols, row1size, row2size, ret_mat)
char *obsfile, *permfile;
int nr, nc, row1size, row2size;
int *cols; /* Ordered column numbers of the datafile to use (some may be repeated) */
int ncols; /* Length of cols vector */
double **ret_mat; /* Must be of size nr by ncols */
{
  int i, j, index=0, *ptri;
  FILE *fpt;
  double *vec, *ptrd, *prow;
  char line[row2size];
  size_t size;

  vec = (double *) malloc(nc*sizeof(double));
  CHECK_MEM(vec);

  size = sizeof(line);

  /* Open the observed p-value file */
  fpt = open_file(obsfile, "r");

  /* Read the first row of column names (not needed) */
  read_header(fpt, row1size);

  /* Read in the p-values */
  fgets(line, size, fpt);
  parse_dline_comma(line, nc, vec);

  fclose(fpt);

  /* Get the -log(p-values) according to the cols vector */
  for (j=0, ptri=cols, ptrd=ret_mat[index]; j<ncols; j++, ptri++, ptrd++) *ptrd = -log(vec[*ptri]);
   
  /* Read the permutation file */
  fpt = open_file(permfile, "r");
  read_header(fpt, row1size);

  /* Read in each line */
  for (i=1; i<=nr; i++) {
    fgets(line, size, fpt);
    parse_dline_comma(line, nc, vec);

    /* Get the p-values according to the cols vector */
    index++;
    prow = ret_mat[index];
    for (j=0, ptri=cols, ptrd=prow; j<ncols; j++, ptri++, ptrd++) *ptrd = vec[*ptri];

    /* Check for missing values */
    if (checkVecBounds(prow, ncols, 0.0, 1.0)) {
      index--;
      continue;
    }

    /* Compute -log(pvalue) */
    for (j=0, ptrd=prow; j<ncols; j++, ptrd++) *ptrd = -log(*ptrd);
  }
  fclose(fpt);
  free(vec);

  /* Return the number of rows of the matrix */
  index++;
  return(index);
  
} /* read_data */

/* Function that R will call */
void ARTP_pathway(obsfile, permfile, nr, nc, cols, ncols, geneStart, geneStop, ngene, rowsize,\
                  inspect_list, llen, iStart, iStop, inspect_path, ilen, r2Files, r2Flag, outfile,\
                  ties_method, method, parms, ret_p_obs, ret_minp, ret_nperm)
char **obsfile, **permfile; /* Data files for observed and permutation p-values */
int *nr, *nc; /* nr is the number of rows of permfile, nc is the number of columns of obsfile and permfile */
int *cols; /* Ordered vector of column numbers from obsfile and permfile to use (some may be repeated) */
int *ncols; /* Length of cols */
int *geneStart, *geneStop; /* The 2 columns of the gene-SNP matrix */
int *ngene;  /* Length of geneStart and geneStop */
int *rowsize; /* 2-element vector for the number of characters of row1and row2 of obsfile (or permfile) */
int *inspect_list, *llen; /* Vector for inspect.gene.list and length of vector */
int *iStart, *iStop; /* Starting and stopping vector positions for each gene for inspect_list. These vectors must be of length ngene. */ 
int *inspect_path, *ilen;  /* Inspect pathway vector */
char **outfile; /* For gene results */
double *ret_p_obs, *ret_minp;  /* Return vectors */
int *ret_nperm, *ties_method;
char **r2Files;
int *r2Flag;
int *method; 
double *parms; 
{
  double **T;
  int nrows, nrp1;

  /* Check for errors */
  fun_ck0(*nc, inspect_list, *llen);

  nrp1 = *nr + 1; 

  /* Allocate memory to store all values */
  T = dblMatrix_alloc(nrp1, *ncols); /* First row is for the observed */

  /* Read the data */
  nrows = read_data(*obsfile, *permfile, *nr, *nc, cols, *ncols, rowsize[0], rowsize[1], T);
  *ret_nperm = nrows - 1;
  if (nrows < 2) return;

  /* For random number generation */
  GetRNGstate();

  /* Call the main function */
  pathway_gene(T, nrows, *ncols, inspect_list, *llen, iStart, iStop, geneStart, geneStop, *ngene,\
                   inspect_path, *ilen, r2Files, *r2Flag, *outfile, ret_p_obs, ret_minp, *ties_method, *method, parms); 

  dblMatrix_free(T, nrp1);

  PutRNGstate();

  return;

} /* END: ARTP_pathway */
 
void R_init_ARTP(DllInfo *dll)
{
    R_registerRoutines(dll, callMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}

