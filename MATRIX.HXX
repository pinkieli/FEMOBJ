//   ********************
//   *** CLASS MATRIX ***
//   ********************
 

#ifndef matrix_hxx

class FloatMatrix ; class IntArray ;


class Matrix
/*
   This abstract class is the superclass of the class that implement
   matrices (FloatMatrix, PolynomialMatrix,...).
 DESCRIPTION :
   A matrix is characterized by its number of rows and columns.
 TASKS :
   Its tasks are defined in the subclasses.
*/
{
   protected:
      int  nRows ;
      int  nColumns ;

   public:
      Matrix ()             { }                         // constructors
      Matrix (int n,int m)  { nRows=n ; nColumns=m ;}
      ~Matrix ()            { }                         // destructor

      void          checkBounds (int,int) ;
      int           giveNumberOfRows ()      { return nRows ;}
      int           giveNumberOfColumns ()   { return nColumns ;}
      int           isSquare ()              { return (! nRows-nColumns) ;}
} ;

#define matrix_hxx
#endif






