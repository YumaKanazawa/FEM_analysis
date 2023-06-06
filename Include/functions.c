#include <stdlib.h>
#include <stdio.h>

#include "functions.h"

/* argument of (x[1],x[2]) in [0,2\pi) or -1 (the origin) */
double arg( double *x )
{
  double eps=1.0e-10;
  double a; 

  if ( x[1]*x[1]+x[2]*x[2] < eps )
    {
      return -1.0;
    }

  if ( x[1] >= fabs(x[2]) )
    {
      a = atan(x[2]/x[1]);
      return (a>=0.0)?a:(a+2*pi);
    }

  if ( x[2] >= fabs(x[1]) )
    {
      a = atan(-x[1]/x[2]);
      return a + 0.5*pi;
    }

  if ( x[1] <= - fabs(x[2]) )
    {
      a = atan(x[2]/x[1]);
      return a+pi;
    }

  if ( x[2] <= - fabs(x[1]) )
    {
      a = atan(-x[1]/x[2]);
      return a + 1.5*pi;
    }
  exit(0);
}


double *dvector(int i, int j) 
{
  double *a;

  if ( (a=(double *)calloc(j-i+1,sizeof(double))) == NULL )
    {
      printf("dvector: memory allocation is failed \n");
      exit(1);
    }  

  return(a-i);
}

void free_dvector(double *a, int i)
{
  free( (void *)(a + i) ); 
}

int *ivector(int i, int j) 
{
  int *a;

  if ( (a=(int *)calloc(j-i+1,sizeof(int))) == NULL )
    {
      printf("ivector: memory allocation is failed \n");
      exit(1);
    }  

  return(a-i);
}

void free_ivector(int *a, int i)
{
  free( (void *)(a + i) ); 
}

double **dmatrix(int nr1, int nr2, int nl1, int nl2)
{
  int i, nrow, ncol; 
  double **a; 

  nrow = nr2 - nr1 + 1 ; /* number of row */
  ncol = nl2 - nl1 + 1 ; /* number of column */

  if ( ( a=(double **)calloc( nrow,sizeof(double *) ) ) == NULL ) 
  {
    printf("dmatrix: memory allocation is failed \n");
    exit(1);
  }
  a = a - nr1; 

  for( i=nr1; i<=nr2; i++) a[i] = (double *)calloc(ncol,sizeof(double)); 
  for( i=nr1; i<=nr2; i++) a[i] = a[i]-nl1;         
  
  return(a);
}

void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2)
{
  int i;

  for ( i = nr1 ; i <= nr2 ; i++) free((void *)(a[i]+nl1));
  free((void *)(a+nr1));
}

int **imatrix(int nr1, int nr2, int nl1, int nl2)
{
  int i, nrow, ncol; 
  int **a; 

  nrow = nr2 - nr1 + 1 ; /* number of row */
  ncol = nl2 - nl1 + 1 ; /* number of column */

  if ( ( a=(int **)calloc( nrow,sizeof(int *) ) ) == NULL ) 
    {
      printf("imatrix: memory allocation is failed \n");
      exit(1);
    }
  a = a - nr1; 

  for( i=nr1; i<=nr2; i++) a[i] = (int *)calloc(ncol,sizeof(int)); 
  for( i=nr1; i<=nr2; i++) a[i] = a[i]-nl1;         
  
  return(a);
}

void free_imatrix(int **a, int nr1, int nr2, int nl1, int nl2)
{
  int i;

  for ( i = nr1 ; i <= nr2 ; i++) free((void *)(a[i]+nl1));
  free((void *)(a+nr1));
}


void matrix_vector_product(double **a, double *b , double *c, int n)
{
  double wk;
  int i, j;

  for ( i = 1; i <= n; i++)
  {
    wk = 0.0;
    for ( j = 1; j <= n; j++ )
    {
      wk += a[i][j]*b[j];
    }
    c[i] = wk;
  }
}

double inner_product( int m, int n, double *a, double *b)
{
  int i;
  double s = 0.0;

  for( i = m; i <= n; i++) s += a[i]*b[i];

  return s ;
}

double vector_norm1( double *a, int m, int n )
{
  int i; 
  double norm = 0.0;
  for ( i = m; i <= n; i++ )
  {
    norm += fabs(a[i]);
  }
  return norm; 
}

double vector_norm2( double *a, int m, int n )
{
  int i; 
  double norm = 0.0;
  for ( i = m; i <= n; i++ )
  {
    norm += a[i]*a[i];
  }
  return sqrt(norm); 
}

void cg(double **a, double *b, double *x, int n, double EPS, int KMAX)
{
  double eps, *r, *p, *tmp, alpha, beta, work; 
  int i, k=0; 

  r = dvector(1,n); /* r[1...n] */
  p = dvector(1,n); /* p[1...n] */
  tmp = dvector(1,n); /* tmp[1...n] */

  matrix_vector_product( a, x, tmp, n );  /* tmp <- A x_0 */

  for( i = 1; i <= n; i++)
  {
    p[i] = b[i] - tmp[i] ; r[i] = p[i];
  }
  /* added by mk */
  eps = vector_norm1(r, 1, n); 
  printf("initial residual=%e\n",eps);
  if ( eps < EPS )  goto OUTPUT;

  do
  { 
    matrix_vector_product( a, p, tmp, n );  /* tmp <- A p_k */
    work = inner_product( 1, n, p, tmp); /* work <- (p,Ap_k) */
    alpha = inner_product( 1, n, p, r ) / work ;

    for( i = 1; i <= n; i++) x[i] = x[i] + alpha*p[i]; 
    for( i = 1; i <= n; i++) r[i] = r[i] - alpha*tmp[i]; 

    k++;  
    eps = vector_norm1(r, 1, n); 
    printf("%d-th residual =%e\n", k, eps);
    if ( eps < EPS )  goto OUTPUT;

    beta = - inner_product( 1, n, r, tmp) / work;
    for( i = 1; i <= n; i++) p[i] = r[i] + beta*p[i];     
 }while( k < KMAX );

 OUTPUT:;

  free_dvector( r, 1 ); free_dvector( p, 1 ); free_dvector( tmp, 1 ); 
  if ( k == KMAX )
  {
    printf("CG did not converge within %d iterations\n", KMAX);
    exit(1);
  }
  else
  {
    printf("CG: %d iterations \n", k); 
  }
}

void matrix_vector_product_crs(double *a,int *ia,int *ja, double *b , double *c, int n)
{
  double wk;
  int i, j;
#pragma omp parallel for private(i,j,wk)
  for ( i = 1; i <= n; i++)
  {
      wk = 0.0;
      for ( j= ia[i]; j < ia[i+1]; ++j )
	{
	  wk += a[j]*b[ja[j]];
	}
      c[i] = wk;
    }
  
}

void cg_crs(double *a,int *ia,int *ja, double *b, double *x, int n, double EPS, int KMAX)
{
  double eps, *r, *p, *tmp, alpha, beta, work; 
  int i, k=0; 

  r = dvector(1,n); /* r[1...n] */
  p = dvector(1,n); /* p[1...n] */
  tmp = dvector(1,n); /* tmp[1...n] */

  matrix_vector_product_crs( a,ia,ja, x, tmp, n );  /* tmp <- A x_0 */

  for( i = 1; i <= n; i++)
  {
    p[i] = b[i] - tmp[i] ; r[i] = p[i];
  }
  /* added by mk */
  eps = vector_norm1(r, 1, n); 
  printf("initial residual=%e\n",eps);
  if ( eps < EPS )  goto OUTPUT;

  do
  { 
    matrix_vector_product_crs( a,ia,ja, p, tmp, n );  /* tmp <- A p_k */
    work = inner_product( 1, n, p, tmp); /* work <- (p,Ap_k) */
    alpha = inner_product( 1, n, p, r ) / work ;

    for( i = 1; i <= n; i++) x[i] = x[i] + alpha*p[i]; 
    for( i = 1; i <= n; i++) r[i] = r[i] - alpha*tmp[i]; 

    k++;  
    eps = vector_norm1(r, 1, n); 
    if(k%100==0) printf("%d-th residual =%e\n", k, eps);
    if ( eps < EPS )  goto OUTPUT;

    beta = - inner_product( 1, n, r, tmp) / work;
    for( i = 1; i <= n; i++) p[i] = r[i] + beta*p[i];     
 }while( k < KMAX );

 OUTPUT:;

  free_dvector( r, 1 ); free_dvector( p, 1 ); free_dvector( tmp, 1 ); 
  if ( k == KMAX )
  {
    printf("CG did not converge within %d iterations\n", KMAX);
    exit(1);
  }
  else
  {
    printf("CG: %d iterations \n", k); 
  }
}

void isort(int *a,int start,int end)
{
  int i,j,tmp;
  for(i=start;i<=end;i++){
    for(j=i+1;j<=end;j++){
      if(a[i]>a[j]){
	tmp=a[i];
	a[i]=a[j];
	a[j]=tmp;
      }
    }
  }
}
