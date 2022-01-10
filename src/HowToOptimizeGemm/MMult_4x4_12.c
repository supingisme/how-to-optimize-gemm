/* Create macros so that the matrices are stored in row-major order */

#define A(i,j) a[ (i)*lda + (j) ]
#define B(i,j) b[ (i)*ldb + (j) ]
#define C(i,j) c[ (i)*ldc + (j) ]

/* Block sizes */
#define mc 256
#define kc 128

#define min( i, j ) ( (i)<(j) ? (i): (j) )

/* Routine for computing C = A * B + C */

void AddDot4x4( int, double *, int, double *, int, double *, int );
void PackMatrixA( int, double *, int, double * );

void MY_MMult( int m, int n, int k, double *a, int lda, 
                                    double *b, int ldb,
                                    double *c, int ldc )
{
  int i, p, pb, ib;

  /* This time, we compute a mc x n block of C by a call to the InnerKernel */

  for ( p=0; p<k; p+=kc ){
    pb = min( k-p, kc );
    for ( i=0; i<m; i+=mc ){
      ib = min( m-i, mc );
      InnerKernel( ib, n, pb, &A( i,p ), lda, &B(p, 0 ), ldb, &C( i,0 ), ldc );
    }
  }
}

void InnerKernel( int m, int n, int k, double *a, int lda, 
                                       double *b, int ldb,
                                       double *c, int ldc )
{
  int i, j;
  double 
    packedA[ m * k ];

  for ( j=0; j<n; j+=4 ){        /* Loop over the columns of C, unrolled by 4 */
    for ( i=0; i<m; i+=4 ){        /* Loop over the rows of C */
      /* Update C( i,j ), C( i,j+1 ), C( i,j+2 ), and C( i,j+3 ) in
	 one routine (four inner products) */
      PackMatrixA( k, &A( i, 0 ), lda, &packedA[ i*k ] );
      AddDot4x4( k, &packedA[ i*k ], 4, &B( 0,j ), ldb, &C( i,j ), ldc );
    }
  }
}

void PackMatrixA( int k, double *a, int lda, double *a_to )
{
  int j;
  
  double
  	*a_0j_pntr = &A(0, 0), *a_1j_pntr = &A(1, 0), 
	*a_2j_pntr = &A(2, 0), *a_3j_pntr = &A(3, 0);

  for( j=0; j<k; j++){  /* loop over columns of A */
    *a_to++ = *a_0j_pntr++;
    *a_to++ = *a_1j_pntr++;
    *a_to++ = *a_2j_pntr++;
    *a_to++ = *a_3j_pntr++;
  }
}



#include <stdio.h>                                                                                                                                                                       
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
 
#include <arm_neon.h>


void AddDot4x4( int k, double *a, int lda,  double *b, int ldb, double *c, int ldc )
{
  /* So, this routine computes a 4x4 block of matrix A

           C( 0, 0 ), C( 0, 1 ), C( 0, 2 ), C( 0, 3 ).  
           C( 1, 0 ), C( 1, 1 ), C( 1, 2 ), C( 1, 3 ).  
           C( 2, 0 ), C( 2, 1 ), C( 2, 2 ), C( 2, 3 ).  
           C( 3, 0 ), C( 3, 1 ), C( 3, 2 ), C( 3, 3 ).  

     Notice that this routine is called with c = C( i, j ) in the
     previous routine, so these are actually the elements 

           C( i  , j ), C( i  , j+1 ), C( i  , j+2 ), C( i  , j+3 ) 
           C( i+1, j ), C( i+1, j+1 ), C( i+1, j+2 ), C( i+1, j+3 ) 
           C( i+2, j ), C( i+2, j+1 ), C( i+2, j+2 ), C( i+2, j+3 ) 
           C( i+3, j ), C( i+3, j+1 ), C( i+3, j+2 ), C( i+3, j+3 ) 
	  
     in the original matrix C 

     A simple rearrangement to prepare for the use of vector registers */

  int p;
  float64x2_t 
    /* hold contributions to
       C( 0, 0 ), C( 0, 1 ), C( 0, 2 ), C( 0, 3 ) 
       C( 1, 0 ), C( 1, 1 ), C( 1, 2 ), C( 1, 3 ) 
       C( 2, 0 ), C( 2, 1 ), C( 2, 2 ), C( 2, 3 ) 
       C( 3, 0 ), C( 3, 1 ), C( 3, 2 ), C( 3, 3 )   */
       c_00_c_01_vreg, c_02_c_03_vreg,
       c_10_c_11_vreg, c_12_c_13_vreg,
       c_20_c_21_vreg, c_22_c_23_vreg,
       c_30_c_31_vreg, c_32_c_33_vreg,
    /* hold 
       A( 0, p ) 
       A( 1, p ) 
       A( 2, p ) 
       A( 3, p ) */
       a_0p_vreg,a_1p_vreg, a_2p_vreg,a_3p_vreg,
       b_p0_b_p1_vreg,
       b_p2_b_p3_vreg;

  double 
    /* Point to the current elements in the four columns of B */
    *a_0p_pntr, *a_1p_pntr, *a_2p_pntr, *a_3p_pntr; 
    
  a_0p_pntr = &A( 0, 0 );
  a_1p_pntr = &A( 1, 0 );
  a_2p_pntr = &A( 2, 0 );
  a_3p_pntr = &A( 3, 0 );

  c_00_c_01_vreg =  vmovq_n_f64(0.f);
  c_02_c_03_vreg =  vmovq_n_f64(0.f);
  c_10_c_11_vreg =  vmovq_n_f64(0.f);
  c_12_c_13_vreg =  vmovq_n_f64(0.f);
  c_20_c_21_vreg =  vmovq_n_f64(0.f);
  c_22_c_23_vreg =  vmovq_n_f64(0.f);
  c_30_c_31_vreg =  vmovq_n_f64(0.f);
  c_32_c_33_vreg =  vmovq_n_f64(0.f);

  for ( p=0; p<k; p++ ){
    b_p0_b_p1_vreg = vld1q_f64(&B(p,0));
    b_p2_b_p3_vreg = vld1q_f64(&B(p,2));
	
    
    a_0p_vreg = vld1q_dup_f64(a_0p_pntr++);
    a_1p_vreg = vld1q_dup_f64(a_1p_pntr++);
    a_2p_vreg = vld1q_dup_f64(a_2p_pntr++);
    a_3p_vreg = vld1q_dup_f64(a_3p_pntr++);


    /* First row and second rows */
    /* c_00_reg += a_0p_reg * b_p0_reg; */
    /* c_01_reg += a_0p_reg * b_p1_reg; */
    /* c_00_c_01_vreg =  vfmaq_laneq_f64(c_00_c_01_vreg, b_p0_b_p1_vreg, a_0p_vreg, 0); */
    c_00_c_01_vreg = vmlaq_f64(c_00_c_01_vreg, b_p0_b_p1_vreg, a_0p_vreg);
    /* c_10_reg += a_1p_reg * b_p0_reg; */
    /* c_11_reg += a_1p_reg * b_p1_reg; */
    /* c_10_c_11_vreg =  vfmaq_laneq_f64(c_10_c_11_vreg, b_p0_b_p1_vreg, a_1p_vreg, 0); */
    c_10_c_11_vreg = vmlaq_f64(c_10_c_11_vreg, b_p0_b_p1_vreg, a_1p_vreg);
    /* c_20_reg += a_2p_reg * b_p0_reg; */
    /* c_21_reg += a_2p_reg * b_p1_reg; */
    /* c_20_c_21_vreg =  vfmaq_laneq_f64(c_20_c_21_vreg, b_p0_b_p1_vreg, a_2p_vreg, 0); */
    c_20_c_21_vreg = vmlaq_f64(c_20_c_21_vreg, b_p0_b_p1_vreg, a_2p_vreg);
    /* c_30_reg += a_3p_reg * b_p0_reg; */
    /* c_31_reg += a_3p_reg * b_p1_reg; */
    /* c_30_c_31_vreg =  vfmaq_laneq_f64(c_30_c_31_vreg, b_p0_b_p1_vreg, a_3p_vreg, 0); */
    c_30_c_31_vreg = vmlaq_f64(c_30_c_31_vreg, b_p0_b_p1_vreg, a_3p_vreg);

    /* c_02_reg += a_0p_reg * b_p2_reg; */
    /* c_03_reg += a_0p_reg * b_p3_reg; */
    /* c_02_c_03_vreg =  vfmaq_laneq_f64(c_02_c_03_vreg, b_p2_b_p3_vreg, a_0p_vreg, 0); */
    c_02_c_03_vreg = vmlaq_f64(c_02_c_03_vreg, b_p2_b_p3_vreg, a_0p_vreg);
    /* c_12_reg += a_1p_reg * b_p2_reg; */
    /* c_13_reg += a_1p_reg * b_p3_reg; */
    /* c_12_c_13_vreg =  vfmaq_laneq_f64(c_12_c_13_vreg, b_p2_b_p3_vreg, a_1p_vreg, 0); */
    c_12_c_13_vreg = vmlaq_f64(c_12_c_13_vreg, b_p2_b_p3_vreg, a_1p_vreg);
    /* c_22_reg += a_2p_reg * b_p2_reg; */
    /* c_23_reg += a_2p_reg * b_p3_reg; */
    /* c_22_c_23_vreg =  vfmaq_laneq_f64(c_22_c_23_vreg, b_p2_b_p3_vreg, a_2p_vreg, 0); */
    c_22_c_23_vreg = vmlaq_f64(c_22_c_23_vreg, b_p2_b_p3_vreg, a_2p_vreg);
    /* c_32_reg += a_3p_reg * b_p2_reg; */
    /* c_33_reg += a_3p_reg * b_p3_reg; */
    /* c_32_c_33_vreg =  vfmaq_laneq_f64(c_32_c_33_vreg, b_p2_b_p3_vreg, a_3p_vreg, 0); */
    c_32_c_33_vreg = vmlaq_f64(c_32_c_33_vreg, b_p2_b_p3_vreg, a_3p_vreg);
  }

  C( 0, 0 ) += c_00_c_01_vreg[0];
  C( 0, 1 ) += c_00_c_01_vreg[1];
  C( 0, 2 ) += c_02_c_03_vreg[0];
  C( 0, 3 ) += c_02_c_03_vreg[1];

  C( 1, 0 ) += c_10_c_11_vreg[0];
  C( 1, 1 ) += c_10_c_11_vreg[1];
  C( 1, 2 ) += c_12_c_13_vreg[0];
  C( 1, 3 ) += c_12_c_13_vreg[1];

  C( 2, 0 ) += c_20_c_21_vreg[0];
  C( 2, 1 ) += c_20_c_21_vreg[1];
  C( 2, 2 ) += c_22_c_23_vreg[0];
  C( 2, 3 ) += c_22_c_23_vreg[1];

  C( 3, 0 ) += c_30_c_31_vreg[0];
  C( 3, 1 ) += c_30_c_31_vreg[1];
  C( 3, 2 ) += c_32_c_33_vreg[0];
  C( 3, 3 ) += c_32_c_33_vreg[1];
}
