/* Create macros so that the matrices are stored in row-major order */

#define A(i,j) a[ (i)*lda + (j) ]
#define B(i,j) b[ (i)*ldb + (j) ]
#define C(i,j) c[ (i)*ldc + (j) ]

/* Routine for computing C = A * B + C */

void AddDot( int, double *, int, double *, double * );

void MY_MMult( int m, int n, int k, double *a, int lda, 
                                    double *b, int ldb,
                                    double *c, int ldc )
{
  int i, j;

  for ( j=0; j<n; j+=1 ){        /* Loop over the columns of C */
    for ( i=0; i<m; i+=1 ){        /* Loop over the rows of C */
      /* Update the C( i,j ) with the inner product of the ith row of A
	 and the jth column of B */
     AddDot( k, &A( i,0 ), ldb, &B( 0,j ), &C( i,j ) );
    }
  }
}


/* Create macro to let X( i ) equal the ith element of x */

#define Y(i) y[ (i)*incy ]

void AddDot( int k, double *x, int incy,  double *y, double *gamma )
{
  /* compute gamma := x' * y + gamma with vectors x and y of length n.

     Here x starts at location x with increment (stride) incx and y starts at location y and has (implicit) stride of 1.
  */
 
  int p;

  for ( p=0; p<k; p++ ){
      /* C(i,j) += A(i,p)*B(p,j) */ 
      /*         => A[ p ] * B( p ); */
      
      /* -> A[p] = A(i,0)[p] =>  x[ p ], x:&A(i, 0) */


      /* -> B(p) = B[p](0,j) -> B[ P*ldb + j ]
       * -> B[ (p)*ldb ](0, j) -> Y(p), y:B(0,j), incy:ldb */

    *gamma += x[ p ] * Y( p ) ;     
  }
}
