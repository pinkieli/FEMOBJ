//   *****************************
//   *** CLASS DIAGONAL MATRIX ***
//   *****************************
 

#include "flotmtrx.hxx"
#include "debug.def"


class DiagonalMatrix : public FloatMatrix
/*
   This class implements a square matrix containing non zero coefficients
   only on the diagonal. Diagonal matrices are typically used as mass
   matrices of finite elements.
 DESCRIPTION :
   A diagonal matrix of size n has nRows=n and nColumns=1. 'values' stores
   only the n diagonal coefficients.
 TASKS :
   Similar to a standard FloatMatrix.
*/
{
   enum  { FALSE , TRUE } ;

   public :
      DiagonalMatrix (int n) : FloatMatrix(n,1)    { }

#     ifdef DEBUG
       double&      at (int i,int j)
			    {this->checkBounds(i,j) ; return values[i-1];}
#     else
       double&      at (int i,int j)               { return values[i-1] ;}
#     endif
      void          checkBounds (int,int) ;
      FloatMatrix*  GiveCopy () ;
      int           isDiagonal ()                  { return TRUE ;}
      void          printYourself () ;
} ;

