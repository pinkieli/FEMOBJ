//   file POLYNOXY.CXX
 
#include "polynoxy.hxx"
#include "flotarry.hxx"
#include <stdlib.h>


PolynomialXY :: PolynomialXY (int n)
   // Constructor. Creates an (X,Y) polynomial of degree n.
{
   degree       = n ;
   coefficients = new FloatArray(++n * ++n / 2) ;
}


double  PolynomialXY :: evaluatedAt (FloatArray* aPoint)
   // Returns the value taken by the receiver at  (x,y) = aPoint .
{
   double answer,x,y,c,cx,cy ;

   x = aPoint -> at(1) ;
   y = aPoint -> at(2) ;

   c      = coefficients -> at(1) ;
   answer = c ;
   if (degree == 0)
      return answer ;

   cx      = coefficients -> at(2) ;
   cy      = coefficients -> at(3) ;
   answer += cx*x + cy*y ;
   if (degree == 1)
      return answer ;

   printf ("(X,Y) polynomials with degree > 1 are not implemented \n") ;
   exit(0) ;
   return 0.;
}


void  PolynomialXY :: printYourself ()
   // Prints the receiver on screen.
{
   printf ("Polynomial(X,Y) of degree %d \n",degree) ;
   coefficients -> printYourself() ;
}








