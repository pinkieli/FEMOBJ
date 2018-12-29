//   file DIAGMTRX.CXX
 
#include "diagmtrx.hxx"


void  DiagonalMatrix :: checkBounds (int i,int j)
   // Checks that the receiver possesses a (i,j) coefficient.
{
   if (i>nRows || i<=0 || j!=i) {
      printf ("error : DiagonalMatrix(%d) has no coeff (%d,%d)\n",nRows,i,j);
      exit(0) ;}
}


FloatMatrix*  DiagonalMatrix :: GiveCopy ()
   // Returns a copy of the receiver. For consistency with method giveCopy()
   // of class FloatMatrix, returns the answer as a FloatMatrix.
{
   DiagonalMatrix *answer ;
   double         *P1,*P2 ;
   int            i ;

   answer = new DiagonalMatrix(nRows) ;
   P1     = values ;
   P2     = answer->values ;
   i      = nRows ;
   while (i--)
      *P2++ = *P1++ ;
   return (FloatMatrix*)answer ;
}


void  DiagonalMatrix :: printYourself ()
   // Prints the receiver on screen.
{
   int i ;

   printf ("Diagonal matrix of size %d\n",nRows) ;
   for (i=0 ; i<nRows ; i++)
      printf ("%f ",values[i]) ;
   printf ("\n") ;
}


