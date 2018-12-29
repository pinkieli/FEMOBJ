//   file FLOTMTRX.CXX

#include "flotmtrx.hxx"
#include "diagmtrx.hxx"
#include "flotarry.hxx"
#include "mathfem.h"

#ifdef DEBUG
double&  FloatMatrix :: at (int i,int j)
   // Returns the coefficient (i,j) of the receiver. Safer but slower than
   // the inline version of method 'at'.
{
   this->checkBounds(i,j) ;
   return values[(j-1)*nRows + i - 1] ;
}
#endif



FloatMatrix*  FloatMatrix :: GiveCopy ()
   // Creates and returns a copy of the receiver.
{
   FloatMatrix *answer ;
   double      *P1,*P2 ;
   int         i ;

   answer = new FloatMatrix(nRows,nColumns) ;
   P1 = answer -> values ;
   P2 = values ;
   i  = nRows * nColumns ;
   while (i--)
      *P1++ = *P2++ ;
   return answer ;
}


double  FloatMatrix :: giveDeterminant ()
   // Returns the determinant of the receiver.
{
#  ifdef DEBUG
      if (! this->isSquare()) {
	 printf ("error : cannot compute determinant of a %d by %d matrix \n"
		 ,nRows,nColumns) ;
	 exit(0) ;}
#  endif

   if (nRows == 1)
      return  values[0] ;
   else if (nRows == 2)
      return  (values[0]*values[3] - values[1]*values[2]) ;
   else {
      printf ("sorry : cannot inverse %d by %d matrices \n",nRows,nColumns) ;
      exit(0) ;}
   return 0.;
}


FloatMatrix*  FloatMatrix :: GiveInverse ()
   // Returns a new matrix, the inverse of the receiver. (implemented only
   // for 1x1 and 2x2 matrices)
{
   FloatMatrix* answer ;
   double       det ;
   double*      p ;

#  ifdef DEBUG
      if (! this->isSquare()) {
	 printf ("error : cannot inverse a %d by %d matrix ! \n",
		  nRows,nColumns);
	 exit(0) ;}
#  endif

   answer = new FloatMatrix(nRows,nRows) ;
   p      = answer->values ;

   if (nRows == 1)
      p[0] = 1. / values[0] ;
   else if (nRows == 2) {
      det  = values[0]*values[3] - values[1]*values[2] ;
      p[0] =  values[3] / det ;
      p[1] = -values[1] / det ;
      p[2] = -values[2] / det ;
      p[3] =  values[0] / det ;}
   else {
      printf ("error : cannot inverse a %d by %d matrix ! \n",nRows,nColumns);
      exit(0) ;}

   return answer ;
}


FloatMatrix*  FloatMatrix :: GiveTransposition ()
   // Returns a new matrix, the transposition of the receiver.
{
   int          i,j ;
   FloatMatrix  *answer ;

   answer = new FloatMatrix(nColumns,nRows) ;
   for (i=1 ; i<=nRows ; i++)
      for (j=1 ; j<=nColumns ; j++)
	 answer->at(j,i) = this->at(i,j) ;
   return answer ;
}



FloatMatrix*  FloatMatrix :: Lumped ()
   // Returns a new diagonal matrix, which is the lumped receiver : all
   // coefficients on a column are concentrated on the diagonal.
{
   DiagonalMatrix *answer ;
   double         sum,*p ;
   int            j ;
   register int   i ;

   answer = new DiagonalMatrix(nRows) ;
   p      = values ;
   for (j=1 ; j<=nColumns ; j++) {
      sum = 0. ;
      i   = nRows ;
      while (i--)
	 sum += *p++ ;
      answer -> at(j,j) = sum ;}

   return answer ;
}


FloatMatrix*  FloatMatrix :: plus (FloatMatrix* aMatrix)
   // Adds aMatrix to the receiver. If the receiver has a null size,
   // adjusts its size to that of aMatrix. Returns the modified receiver.
{
   register int i ;
   int      n,m ;
   double   *P1,*P2 ;

   if (aMatrix -> isDiagonal())
      return  this->plusDiagonalMatrix((DiagonalMatrix*)aMatrix) ;

   n = aMatrix -> nRows ;
   m = aMatrix -> nColumns ;
   if (nRows*nColumns == 0) {
      if (values)
	 freeDouble (values) ;
      nRows    = n ;
      nColumns = m ;
      i        = n * m ;
      values   = allocDouble(i) ;
      P1       = values ;
      P2       = aMatrix->values ;
      while (i--)
	 *P1++ = *P2++ ;}
   else {
#     ifdef DEBUG
	 if (n-nRows || m-nColumns) {
	     printf ("dimensions mismatch : r1,c1,r2,c2 : %d %d %d %d\n",
		      nRows,n,nColumns,m) ;
	     exit(0) ;}
#     endif

      P1 = values ;
      P2 = aMatrix->values ;
      i  = n * m ;
      while (i--)
	 *P1++ += *P2++ ;}

   return this ;
}


FloatMatrix*  FloatMatrix :: plusDiagonalMatrix (DiagonalMatrix* aMatrix)
{
   register int i ;
   int      n ;

   n = aMatrix -> giveNumberOfRows() ;
   if (nRows*nColumns == 0) {
      if (values)
	 freeDouble (values) ;
      nRows    = n ;
      nColumns = n ;
      values   = allocDouble(n*n) ;}

#  ifdef DEBUG
      if (n-nRows) {
	 printf ("dimensions mismatch in FloatMatrix+DiagonalMatrix\n") ;
	 exit(0) ;}
#  endif

   for (i=1 ; i<=nRows ; i++)
      this->at(i,i) += aMatrix->at(i,i) ;

   return this ;
}


void  FloatMatrix :: plusProduct (FloatMatrix* a, FloatMatrix* b, double dV)
   // Adds to the receiver the product  a(transposed).b dV .
   // If the receiver has a null size, it is expanded.
   // This method assumes that both the receiver and the product above are
   // symmetric matrices, and therefore computes only the upper half of the
   // receiver ; the lower half is not modified. Other advantage : it does
   // not compute the transposition of matrix a.
{
   int    i,j,p ;
   double coeff ;

   double *P1,*P2 ;

   if (nRows*nColumns == 0) {
      nRows  = nColumns = a->nColumns ;
      values = allocDouble(nRows*nColumns) ;}

   p = a->nRows ;
   for (i=1 ; i<=nRows ; i++)
      for (j=i ; j<=nColumns ; j++) {
	 P1 = &(a->at(1,i)) ;
	 P2 = &(b->at(1,j)) ;
	 //coeff = dotProduct(P1,P2,p) ;
	 coeff = 0.;
	 while (p--)
        coeff += *P1++ * *P2++ ;
     //
	 this->at(i,j) += coeff * dV ;}
}


void  FloatMatrix :: printYourself ()
   // Prints the receiver on screen.
{
   int i,j ;

   printf ("FloatMatrix with dimensions : %d %d\n",
	    nRows,nColumns) ;
   if (nRows<=10 && nColumns<=8)
      for (i=1 ; i<=nRows ; ++i) {
	 for (j=1 ; j<=nColumns && j<10 ; ++j)
	    printf ("%f  ",this->at(i,j)) ;
	 printf ("\n") ; }
   else
      printf ("   large matrix : coefficients not printed \n") ;
}


FloatMatrix*  FloatMatrix :: rotatedWith (FloatMatrix* r)
   // Returns the receiver 'a' rotated according the change-of-base matrix r.
   // The method performs the operation  a = r(transp) * a * r .
{
   register int  i ;
   double        *p1,*p2 ;
   FloatMatrix   *rt,*rta,*rtar ;

   rt   = r -> GiveTransposition() ;          //  r(transp)
   rta  = rt -> Times(this) ;                 //  r(transp) . a
   rtar = rta -> Times(r) ;                   //  r(transp) . a . r

   p1 = values ;
   p2 = rtar -> values ;
   i  = nRows*nColumns ;
   while (i--)
      *p1++ = *p2++ ;                         // copy rtar into the receiver

   delete rt ;
   delete rta ;
   delete rtar ;
   return this ;
}


FloatMatrix*  FloatMatrix :: symmetrized ()
   // Initializes the lower half of the receiver to the upper half.
{
   int i,j ;

#  ifdef DEBUG
      if (nRows != nColumns) {
	 printf ("error : cannot symmetrize a non-square matrix \n") ;
	 exit(0) ;}
#   endif

   for (i=2 ; i<=nRows ; i++)
      for (j=1 ; j<i ; j++)
	 this->at(i,j) = this->at(j,i) ;

   return this ;
}


FloatMatrix*  FloatMatrix :: times (double factor)
   // Multiplies every coefficient of the receiver by factor. Answers the
   // modified receiver.
{
   register int i ;
   double*      p ;

   p = values ;
   i = nRows * nColumns ;
   while (i--)
      *p++ *= factor ;
   return this ;
}


FloatArray*  FloatMatrix :: Times (FloatArray* anArray)
   // Returns the product of the receiver and anArray.
{
   int         i,j ;
   double      sum ;
   FloatArray* answer ;

#  ifdef DEBUG
      if (nColumns - anArray->giveSize()) {
	 printf ("error in product A*x : dimension mismatch \n") ;
	 this    -> printYourself() ;
	 anArray -> printYourself() ;
	 exit(0) ;}
#  endif

   answer = new FloatArray(nRows) ;
   for (i=1 ; i<=nRows ; i++) {
      sum = 0. ;
      for (j=1 ; j<=nColumns ; j++)
	 sum += this->at(i,j) * anArray->at(j) ;
      answer->at(i) = sum ;}
   return answer ;
}


FloatMatrix*  FloatMatrix :: Times (FloatMatrix* aMatrix)
   // Returns the product of the receiver and aMatrix. Easier to use than
   // operator * .
{
   int          i,j,k,p ;
   double       coeff ;
   FloatMatrix* answer ;

#  ifdef DEBUG
      if (nColumns != aMatrix->nRows) {
	 printf ("error in product A*B : dimensions do not match \n") ;
	 this    -> printYourself() ;
	 aMatrix -> printYourself() ;
	 exit(0) ;}
#  endif

   p      = aMatrix -> nColumns ;
   answer = new FloatMatrix(nRows,p) ;
   for (i=1 ; i<=nRows ; i++)
      for (j=1 ; j<=p ; j++) {
         coeff = 0. ;
	 for (k=1 ; k<=nColumns ; k++)
	    coeff += this->at(i,k) * aMatrix->at(k,j) ;
	 answer->at(i,j) = coeff ;}
   return answer ;
}



