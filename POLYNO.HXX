//   ************************
//   *** CLASS POLYNOMIAL ***
//   ************************


#ifndef polyno_hxx

#include "flotarry.hxx"
#include <stdio.h>


class Polynomial
/*
   This abstract class is the superclass of the classes that implement
   polynomials : polynomials with one variable (X), with two variables
   (X,Y), etc.
 DESCRIPTION :
   The degree of a polynomial is stored in attribute 'degree ; for example,
   P(X,Y) = 2 + 4 X + 5 XY has degree 2.
   The coefficients of the polynomial are stored in 'coefficients' ; their
   sequence depends on the the polynomial's type (PolynomialXY, etc).
 TASKS :
   The task of a polynomial is to return its value at a given point (method
   'evaluatedAt').
*/
{
   protected :
      int          degree ;
      FloatArray*  coefficients ;

   public :
      Polynomial ()                       { }
      virtual ~Polynomial ()              { delete coefficients;}

      double&         at (int i)          { return coefficients -> at(i) ;}
      virtual double  evaluatedAt (FloatArray*)  = 0 ; // 0 => NUll pinkie.ljz
      virtual void    printYourself ()           = 0 ; // 0 => NUll pinkie.ljz
} ;

#define polyno_hxx
#endif








