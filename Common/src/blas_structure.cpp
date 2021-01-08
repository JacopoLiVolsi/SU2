/*!
 * \file blas_structure.cpp
 * \brief Implementation of the functions that either simulate BLAS functionality
          or interface to an actual BLAS implementation.
 * \author E. van der Weide
 * \version 6.0.1 "Cardinal"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/blas_structure.hpp"
#include <cstring>

/* MKL or BLAS, if supported. */
#if (defined (HAVE_MKL) || defined(HAVE_BLAS)) && !(defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))

/* Function prototypes for the BLAS routines used. */
extern "C" void dgemm_(char*, char*, const int*, const int*, const int*,
                       const passivedouble*, const passivedouble*, const int*,
                       const passivedouble*, const int*,
                       const passivedouble*, passivedouble*, const int*);

extern "C" void dgemv_(char*, const int*, const int*, const passivedouble*,
                       const passivedouble*, const int*, const passivedouble*,
                       const int*, const passivedouble*, passivedouble*, const int*);

extern "C" void dgetrf_(const int*, const int*, passivedouble*, const int*,
                              int*,       int*);

extern "C" void dgetrs_(char*, const int*, const int*, const passivedouble*,
                        const int*, int*, passivedouble*, const int*, int*);

extern "C" void dgels_(char*, const int*, const int*, const int*, passivedouble*, const int*, 
		       passivedouble*, const int*, passivedouble*, const int*, const int*);

#endif

/* Constructor. Initialize the const member variables, if needed. */
CBlasStructure::CBlasStructure(void)
#if !(defined(HAVE_LIBXSMM) || defined(HAVE_BLAS) || defined(HAVE_MKL)) || (defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))
  : mc (256), kc (128), nc (128) 
#endif
{}

/* Destructor. Nothing to be done. */
CBlasStructure::~CBlasStructure(void) {}

/* Dense matrix multiplication, gemm functionality. */
void CBlasStructure::gemm(const int M,        const int N,        const int K,
                          const su2double *A, const su2double *B, su2double *C,
                          CConfig *config) {

  /* Initialize the variable for the timing, if profiling is active. */
#ifdef PROFILE
  double timeGemm;
  if( config ) config->GEMM_Tick(&timeGemm);
#endif

#if (defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)) || !(defined(HAVE_LIBXSMM) || defined(HAVE_MKL) || defined(HAVE_BLAS))
  /* Native implementation of the matrix product. This optimized implementation
     assumes that the matrices are in column major order. This can be
     accomplished by swapping N and M and A and B. This implementation is based
     on https://github.com/flame/how-to-optimize-gemm. */
  gemm_imp(N, M, K, B, A, C);

#else
#ifdef HAVE_LIBXSMM

  /* The gemm function of libxsmm is used to carry out the multiplication.
     Note that libxsmm_gemm expects the matrices in column major order. That's
     why the in the calling sequence A and B and M and N are reversed. */
  su2double alpha = 1.0;
  su2double beta  = 0.0;
  char trans = 'N';

  libxsmm_dgemm(&trans, &trans, &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);

#else // MKL and BLAS

  /* The standard blas routine dgemm is used for the multiplication.
     Call dgemm without transposing the matrices. In that case dgemm expects
     the matrices in column major order, see the comments for libxsmm. */
  su2double alpha = 1.0;
  su2double beta  = 0.0;
  char trans = 'N';

  dgemm_(&trans, &trans, &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);

#endif
#endif

  /* Store the profiling information, if needed. */
#ifdef PROFILE
  if( config ) config->GEMM_Tock(timeGemm, M, N, K);
#endif
}

/* Dense matrix vector multiplication, gemv functionality. */
void CBlasStructure::gemv(const int M,        const int N,   const su2double *A,
                          const su2double *x, su2double *y) {

#if (defined (HAVE_BLAS) || defined(HAVE_MKL)) && !(defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))

  /* The standard blas routine dgemv is used for the multiplication.
     Note that dgemv expects the matrices in column major order, while
     A is in row major order. This can be solved by using the transpose
     and switching M and N. */
     su2double alpha = 1.0;
     su2double beta  = 0.0;
     int       inc   = 1;
     char trans = 'T';

     dgemv_(&trans, &N, &M, &alpha, A, &N, x, &inc, &beta, y, &inc);

#else

  /* Native implementation of the matix vector product.
     Initialize the elements of y to zero. */
  memset(y, 0, M*sizeof(su2double));  

  /* Carry out the matrix vector product. */
  for(int k=0; k<M; ++k) {
    const su2double *AA = A + k*N;
    for(int l=0; l<N; ++l)
      y[k] += AA[l]*x[l];
  }

#endif
}

#if !(defined(HAVE_LIBXSMM) || defined(HAVE_BLAS) || defined(HAVE_MKL)) || (defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))

/* Macros for accessing submatrices of a matmul using the leading dimension. */
#define A(i, j) a[(j)*lda + (i)]
#define B(i, j) b[(j)*ldb + (i)]
#define C(i, j) c[(j)*ldc + (i)]

/* Function, which perform the implementation of the gemm functionality.  */
void CBlasStructure::gemm_imp(const int m,        const int n,        const int k,
                              const su2double *a, const su2double *b, su2double *c) {

  /* Initialize the elements of c to zero. */
  memset(c, 0, m*n*sizeof(su2double));

  /* Set the leading dimensions of the three matrices. */
  const int lda = m;
  const int ldb = k;
  const int ldc = m;

  /* The full matrix multiplication is split in several blocks.
     Loop over these blocks. */
  for(int p=0; p<k; p+=kc) {
    int pb = min(k-p, kc);
    for(int j=0; j<n; j+=nc) {
      int jb = min(n-j, nc);
      for(int i=0; i<m; i+=mc) {
        int ib = min(m-i, mc);

        /* Carry out the multiplication for this block. */
        gemm_inner(ib, jb, pb, &A(i, p), lda, &B(p, j), ldb, &C(i, j), ldc);
      }
    }
  }
}

/* Compute a portion of the c matrix one block at a time.
   Handle ragged edges with calls to a slow but general function. */
void CBlasStructure::gemm_inner(int m, int n, int k, const su2double *a, int lda,
                                const su2double *b, int ldb, su2double *c, int ldc) {

  /* Carry out the multiplication for this block. At the
     moment simply a call to gemm_arbitrary. */
  gemm_arbitrary(m, n, k, a, lda, b, ldb, c, ldc);
}

/* Naive gemm implementation to handle arbitrary sized matrices. */
void CBlasStructure::gemm_arbitrary(int m, int n, int k, const su2double *a, int lda,
                                    const su2double *b, int ldb, su2double *c, int ldc) {

  /* The order of these loops is tuned for column-major matrices. */
  for (int p = 0; p < k; p++) {
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        C(i, j) += A(i, p) * B(p, j);
      }
    }
  }
}

#undef C
#undef B
#undef A

#endif


void CBlasStructure::dgetrs(const int N, const int NRHS, su2double *A, su2double *x){

 const int lda = N;
 const int ldb = NRHS;
 int ii, jj, kk;
 int p[N];
 int Rank = N;
 for(ii = 0; ii < N; ii++)
  p[ii] = ii;
 su2double pivot = 0.0, maxpivot = 0.0;
 int imaxpivot = 0;
 bool null_pivot = false;

 su2double det = 1.0;
 su2double det_toll = 1E-2, toll = 1E-4;
/*- Compute the infinity norm of A -*/
 su2double b_norm = Norm2(NRHS, x);
 su2double A_norm_inf = 0.0;
 for(jj = 0; jj < N; jj++)
  A_norm_inf += abs(A[jj*N]);
 for(ii = 1; ii < N; ii++){
  su2double aaa = 0.0;
  for(jj = 0; jj < N; jj++)
   aaa += abs(A[ii + jj*N]);
  A_norm_inf = max(A_norm_inf,aaa);
 }
/*- Copy to keep the non-factorized A -*/
 su2double *A2 = new su2double[N*N];  
 for(ii = 0; ii < N*N; ii++) A2[ii] = A[ii];

#if (defined (HAVE_MKL) || defined(HAVE_BLAS)) && !(defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))
/* MKL or BLAS, if supported. */
 int info, work;
 char NoTrans = 'N';
 int info_solve;
 su2double *work_vec;
/*--- Factorize matrix A ---*/ 
 dgetrf_(&lda, &lda, A, &lda, &piv, &info);

 for(ii = 0; ii < N; ii++)
  det *= A[ii+N*ii];

if(abs(det)/A_norm_inf > det_toll)
/*--- Solve linear system L*U*x = b ---*/ 
 dgetrs_(&NoTrans, &lda, &ldb, A, &lda, &piv, x, &ldb, &info_solve);
else
  dgels_(&Norans,&ldb,&lda,&lda,A2,&ldb,x,&lda,work_vec,&work,&info);

#else
 
/*--- LU factorization with pivoting ---*/ 
 for(ii = 0; ii < N-1; ++ii){ 
  maxpivot = A[p[ii] + ii*lda];
  imaxpivot = ii;
  for(kk = ii+1; kk < N; ++kk){ 
   if(A[p[kk]+ ii*lda] > maxpivot){
    maxpivot = A[p[kk]+ ii*lda];
    imaxpivot = kk;
   }
  }

  if(imaxpivot != ii)
   swap (p[ii], p[imaxpivot]);      
       
  pivot = A[p[ii] + ii*lda]; 

/*--- Check for singular matrices  ---*/
  if(abs(pivot)/b_norm < toll && Rank > 1){
   Rank--;
   null_pivot = true;
  } 

  for(jj = ii+1; jj < N; jj++){
   if(!null_pivot)
    A[p[jj] + ii*lda] = A[p[jj] + ii*lda] / pivot;       
   for(kk = ii+1; kk < N; kk++)
    A[p[jj] + kk*lda] += -1.0*A[p[ii]+ kk*lda]*A[p[jj] + ii*lda];
  }

   null_pivot = false;        

 } // end for ii

 for(ii = 0; ii < N; ii++)
  det *= A[ii+N*ii];

//cout<<"Determinant is "<< det;
if(abs(det)/A_norm_inf > det_toll){

cout<<", LU method used."<<endl;
/*--- Solve linear system L*U*x = b ---*/ 

 su2double  f = 0.0;
 su2double *y = new su2double[NRHS];
 /*- Forward substitution -*/
 for(ii = 0; ii < N; ii++){
  f = x[p[ii]];
  for (kk = 0; kk < ii; kk++)
   f -= A[p[ii] + kk*lda] * y[p[kk]];
  y[p[ii]] = f;
 }

 /*- Backward substitution -*/
 for(jj = 1; jj <= N; jj++){
  ii = N - jj;
  f = y[p[ii]];
  for (kk = ii+1; kk < N; kk++)
   f -= A[p[ii] + kk*lda]*x[p[kk]];
  if(abs(A[p[ii] + ii*lda]) > toll)
   x[p[ii]] = f / A[p[ii] + ii*lda]; 
  else
   x[p[ii]] =  0.0;
 }

 /*- Backward substitution -
  for(ii = N-1; ii != -1u; ii--){
   f = y[p[ii]];
   for(jj = N-1; jj > ii; jj--)
    f -= A[p[ii] + jj*lda]*x[p[jj]];
  if(abs(A[p[ii] + ii*lda]) > toll)
   x[p[ii]] = f / A[p[ii] + ii*lda]; 
  else
   x[p[ii]] =  0.0;
  }*/
    
 delete [] y;

} else {

cout<<", Rank is "<<Rank<<" N is "<<N<<", LS method used."<<endl;
  su2double *A_rect, *b;
  A_rect = new su2double[Rank*N];
  b      = new su2double[Rank]; 
  for(ii = 0; ii < Rank; ii++)
   for(jj = 0; jj < N; jj++)
    A_rect[ii + jj*Rank] = A2[p[ii] + jj*N];
  for(jj = 0; jj < Rank; jj++)
   b[jj] = x[p[jj]];

  dgels(Rank, N, A_rect, b, x);

  delete [] b;
  delete [] A_rect;

 } // end else

#endif

 delete [] A2;

 return;
} // end function
 

void CBlasStructure::dgels(const int Mrow, const int Ncol, su2double* A, su2double* b, su2double* x){

 int ii, jj, kk;
 
/*--- QR factorization of A^{T} \in |R^{Ncol x Mrow}, with Ncol > Mrow ---*/
 su2double *Q = new su2double[Ncol*Ncol]; // Q \in |R^{Ncol x Ncol}
 su2double *R = new su2double[Ncol*Mrow]; // R \in |R^{Ncol x Mrow}

 for(ii = 0; ii < Ncol*Ncol; ii++) Q[ii] = 0.0;
 for(ii = 0; ii < Ncol*Mrow; ii++) R[ii] = 0.0;

 su2double *a_i = new su2double[Ncol];
 su2double *q_i = new su2double[Ncol];

 su2double *y   = new su2double[Ncol];
 for(ii = 0; ii < Ncol; ii++)
  y[ii] = 0.0;

 su2double toll = 1E-5;
 su2double A_norm_inf = 0.0;
 for(jj = 0; jj < Ncol; jj++)
  A_norm_inf += abs(A[jj*Mrow]);
 for(ii = 1; ii < Ncol; ii++){
  su2double aaa = 0.0;
  for(jj = 0; jj < Mrow; jj++)
   aaa += abs(A[ii + jj*Mrow]);
  A_norm_inf = max(A_norm_inf,aaa);
 }

 su2double aux_norm;
 su2double aux_prod; 

/*-- Compute matrix Q -> Modified Gram-Schmidt --*/

// Initializing the first column of Q 
 for(ii = 0; ii < Ncol; ii++)
  a_i[ii] = A[ii*Mrow]; 
 aux_norm = Norm2(Ncol, a_i); //if(aux_norm == 0.0)cout<<"NORMA NULLA."<<endl;
 for(ii = 0; ii < Ncol; ii++)
  Q[ii] = a_i[ii]/aux_norm;  

// Computing the remaining column of Q
 for(jj = 1; jj < Mrow; jj++){

  for(ii = 0; ii < Ncol; ii++)
   a_i[ii] = A[jj + ii*Mrow];

  for(kk = 0; kk < jj; kk++){

   for(ii = 0; ii < Ncol; ii++)
    q_i[ii] = Q[ii + Ncol*kk];   

   aux_prod = Scalar(a_i,q_i,Ncol);
   for(ii = 0; ii < Ncol; ii++)
    a_i[ii] -= aux_prod*q_i[ii];

  } // end kk

  aux_norm = Norm2(Ncol, a_i); //if(aux_norm == 0.0)cout<<"NORMA NULLA."<<endl;
  if(aux_norm != 0.0){
   for(ii = 0; ii < Ncol; ii++)
     Q[ii + Ncol*jj] = a_i[ii]/aux_norm;
  }

 } // end jj

//cout<<"Orthogonal Q"<<endl;for(ii = 0; ii < Ncol; ii++){cout<<"[ ";for(jj = 0; jj < Ncol; jj++)cout<<Q[ii+jj*Ncol]<<" ";cout<<"]"<<endl;}
//cout << "Determinat of Q = " << Determinant(Q, Ncol) << endl;
/*-- Compute matrix R = Q^{T}*A^{T} --*/
// gemm_imp(Ncol, Mrow, Mrow, Q, A, R);

 for(ii = 0; ii < Mrow; ii++){
  for(jj = 0; jj < Mrow; jj++){
   aux_prod = 0.0;
   for(kk = 0; kk < Ncol; kk++){
    aux_prod += Q[kk + ii*Ncol]*A[jj + Mrow*kk]; // Q[ii + kk*Ncol] -> sbagliato?
   }
   R[ii + Ncol*jj] = aux_prod;
  }
 }
//cout<<"Is R triang sup?"<<endl;for(ii = 0; ii < Ncol; ii++){cout<<"[ ";for(jj = 0; jj < Mrow; jj++)cout<<R[ii+jj*Ncol]<<" ";cout<<"]"<<endl;}

/*--- System R^{T}y = b solution. ---*/
 su2double f;
 /*- Forward substitution -*/
 for(ii = 0; ii < Mrow; ++ii){
  f = b[ii];
  for(kk = 0; kk < ii; ++kk){
    f -= R[kk + ii*Ncol] * y[kk];
  }
  if(abs(R[ii + ii*Ncol])/A_norm_inf > toll) 
   y[ii] = f/R[ii + ii*Ncol];
  else
   y[ii] = 0.0;
 }

/*--- System Q^{T}x = y solution. ---*/
 for(ii = 0; ii < Ncol; ii++){
  x[ii] = 0.0;
  for(jj = 0; jj < Ncol; jj++)
   x[ii] += Q[jj + Ncol*ii]*y[jj];
 }

 delete [] Q;   delete [] R;
 delete [] a_i; delete [] q_i;
 delete [] y;

 return;
}


void CBlasStructure::Iterative_Solve(const int N, const int NRHS, su2double *A, su2double *x){

 int it, ii, jj;
 it = 0;
 const int max_it = 99; 
 const su2double toll = 1.0E-6;
 su2double  res_norm = 1, res_0_norm;
 su2double *x_0 = new su2double[N];
 su2double *res = new su2double[N]; 
 su2double *z   = new su2double[N]; 
 su2double *aux = new su2double[N]; 
 for(ii = 0; ii < N; ii++){
  x_0[ii] = 0.0;
  res[ii] = 0.0;
  z[ii]   = 0.0;
  aux[ii] = 0.0;
 }

 su2double *P_Jacobi_inv = new su2double[N];
 su2double *AT = NULL; 

// Hp: A not symmetric

 AT = new su2double[N*N];
 for(ii = 0; ii < N; ii++){
  su2double ppp = 0.0;
  for(jj = 0; jj < N; jj++)
   ppp += A[jj*N + ii]*A[jj*N + ii];
  P_Jacobi_inv[ii] = 1/sqrt(ppp);
 }
 for(ii = 0; ii < N; ii++){
  for(jj = 0; jj < N; jj++)
   AT[jj + ii*N] = A[ii + jj*N];
 }
 gemv(N, N, AT, x_0, res);

//
/* Hp: A symmetric
 for(ii = 0; ii < N; ii++){
  P_Jacobi_inv[ii] = 1/A[ii + N*ii];
 }
 gemv(N, N, A, x_0, res);
*/

// Metodo Jacobi
 for(ii = 0; ii < N; ii++)
  res[ii] = x[ii] - res[ii];
 res_0_norm = Norm2(N, res);
 if(res_0_norm == 0.0)
  res_0_norm = 1.0;

 while( (it <= max_it) && (res_norm > toll) ){

  for(ii = 0; ii < N; ii++){
   z[ii]    = P_Jacobi_inv[ii]*res[ii];
   x_0[ii] += z[ii];
  }

  gemv(N, N, A, z, aux);
  for(ii = 0; ii < N; ii++)
   res[ii] -= aux[ii];
 
// Normalized residual criteria
  res_norm  = Norm2(N, res);
  res_norm /= res_0_norm;
  it++;
 }

 for(ii = 0; ii < N; ii++)
  x[ii] = x_0[ii];

 if(x_0 != NULL) delete [] x_0;
 if(aux != NULL) delete [] aux;
 if(z   != NULL) delete [] z;
 if(res != NULL) delete [] res;
 if(AT  != NULL) delete [] AT;
 if(P_Jacobi_inv != NULL) delete [] P_Jacobi_inv;

return;
} // end function


void CBlasStructure::Conjugate_Gradient(const int N, const int NRHS, const su2double *A, su2double *x){

 int it, ii, jj;
 it = 0;
 const int max_it = 99; 
 const su2double toll = 1.0E-6;
 su2double  res_norm = 1, res_norm_0;
 su2double *xk  = new su2double[N];
 su2double *pk  = new su2double[N]; 
 su2double *rk  = new su2double[N]; 
 su2double *Apk = new su2double[N]; 
 su2double *Axk = new su2double[N]; 
 su2double *Ark = new su2double[N]; 
 su2double  pkApk, alpha, beta;

/* A need to be row-major order for gemv(), but A is symmetric */

 for(ii = 0; ii < N; ii++){
  xk[ii] = 0.0;
  rk[ii] = x[ii];
  pk[ii] = rk[ii];
 }
 res_norm_0 = Norm2(N, rk);
 if(res_norm_0 == 0.0)
  res_norm_0 = 1.0;

 while( (it <= max_it) && (res_norm > toll) ){

// Compute A*pk and the product (pk,A*pk)
  for(ii = 0; ii < N; ii++){
   Apk[ii] = 0.0;
   for(jj = 0; jj < N; jj++)
    Apk[ii] += A[ii + jj*N]*pk[jj];
  }
//  gemv(N, N, A, pk, Apk);
  pkApk = Scalar(pk,Apk,N);

// Compute alpha
  alpha = 0.0;
  if(pkApk != 0.0)
   alpha = Scalar(rk,pk,N)/pkApk; 

// Update xk
  for(ii = 0; ii < N; ii++)
   xk[ii] = xk[ii] + alpha*pk[ii];

// Update rk
  for(ii = 0; ii < N; ii++){
   Axk[ii] = 0.0;
   for(jj = 0; jj < N; jj++)
    Axk[ii] += A[ii + jj*N]*xk[jj];
  }
//  gemv(N, N, A, xk, Axk);
  for(ii = 0; ii < N; ii++)
   rk[ii] = rk[ii] - alpha*Axk[ii];

// Compute beta
  for(ii = 0; ii < N; ii++){
   Ark[ii] = 0.0;
   for(jj = 0; jj < N; jj++)
    Ark[ii] += A[ii + jj*N]*rk[jj];
  }
//  gemv(N, N, A, rk, Ark);
  beta = 0.0;
  if(pkApk != 0.0){
   beta = Scalar(pk,Ark,N)/pkApk;
//   beta = Scalar(rk,Apk,N)/pkApk;
  }

// Update pk
  for(ii = 0; ii < N; ii++)
   pk[ii] = rk[ii] - beta*pk[ii];

// Convergence check with normalized residual criteria
  res_norm  = Norm2(N, rk);
  res_norm /= res_norm_0;
  it++;

 } // end while

/*--- Update for output ---*/
 for(ii = 0; ii < N; ii++)
  x[ii] = xk[ii];

 if(xk  != NULL) delete [] xk;
 if(rk  != NULL) delete [] rk;
 if(pk  != NULL) delete [] pk;
 if(Apk != NULL) delete [] Apk;
 if(Axk != NULL) delete [] Axk;
 if(Ark != NULL) delete [] Ark;

} // end function




su2double CBlasStructure::Norm2(const int N, su2double* x){

 su2double sum = 0.0;

 for(unsigned int ii = 0; ii < N; ii++)
  sum += x[ii]*x[ii];

 return sqrt(sum);
}

su2double CBlasStructure::Determinant(su2double* A, int N){

 su2double det = 0.0;
 if(N==1){
  det = A[0];
 }
 else if(N==2){
  det = A[0]*A[3] - A[1]*A[2];
 }
 else{
  for(int row = 0; row < N; row++){
   su2double sub_A[N*N];
    for (int ii = 0; ii < N-1; ii++) {
     for (int jj = 0; jj < N-1; jj++) {
      int sub_ii = (ii < row ? ii : ii+1);
      int sub_jj = jj+1;
      sub_A[jj*N+ii] = A[sub_jj*N+sub_ii];
     }
    }
   if(row%2 == 0)
    det += A[row]*Determinant(sub_A, N-1);
   else
    det -= A[row]*Determinant(sub_A, N-1);
  }
 }

 return det;
} // end function


su2double CBlasStructure::Scalar(su2double *a, su2double *b, int N){
 su2double res = 0.0;
 for(int ii = 0; ii < N; ii++)
  res += a[ii]*b[ii];
 return res;
} // end function


