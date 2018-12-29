//   file POLYMTRX.CXX

#include "polymtrx.hxx"
#include "polyno.hxx"
#include "flotmtrx.hxx"


PolynomialMatrix :: PolynomialMatrix (int n, int m)
		  : Matrix (n,m)
   // Constructor. Creates a n-by-m matrix of polynomials.
{
   values = new Polynomial* [n*m] ;
   this -> initialized () ;
}


PolynomialMatrix :: ~PolynomialMatrix ()
   // Destructor.
   // Warning : if the receiver contains several times the same polynomial,
   //           or if several matrices contain the same polynomial, multiple
   //           attempts to destroy the polynomial may be desastrous !
   //           a possible solution may be for matrices to store only copies
   //           of the polynomials.
{
   int i,size ;

   i = size = nRows*nColumns ;
   while (i--)
      delete values[i] ;
   delete [] values ;
}


#ifdef DEBUG
Polynomial*&  PolynomialMatrix :: at (int i, int j)
   // Returns a polynomial, the coefficient (i,j) of the receiver. Slow but
   // safe.
{
   this -> checkBounds(i,j) ;
   return  values[(j-1)*nRows + i -1] ;
}
#endif


FloatMatrix*  PolynomialMatrix :: EvaluatedAt (FloatArray* aPoint)
   // Returns a matrix of numbers, which is the receiver where all polynomial
   // coefficients are evaluated at the point aPoint.
{
   FloatMatrix* answer ;
   int          i,j ;

   answer = new FloatMatrix(nRows,nColumns) ;
   for (i=1 ; i<=nRows ; i++)
      for (j=1 ; j<=nColumns ; j++)
	 answer -> at(i,j) = this -> at(i,j) -> evaluatedAt(aPoint) ;
   return answer ;
}


PolynomialMatrix*  PolynomialMatrix :: initialized ()
   // Sets all values of the receiver to nil. Returns the receiver.
{
   register int  i ;
   Polynomial**  p ;

   p = values ;
   i = nRows * nColumns ;
   while (i--)
      *p++ = NULL ;
   return this ;
}


#ifdef DEBUG
void  PolynomialMatrix :: printYourself ()
   // Prints the receiver on screen.
{
   int         i,j ;
   Polynomial* p ;

   printf ("PolynomialMatrix of size [%d,%d] \n",nRows,nColumns) ;
   for (i=1 ; i<=nRows ; i++)
      for (j=1 ; j<=nColumns ; j++) {
	 printf ("  (%d,%d) :  ",i,j) ;
	 p = this -> at(i,j) ;
	 if (p)
	   p -> printYourself() ;
	 else
	   printf ("NULL \n") ;}
}
#endif


