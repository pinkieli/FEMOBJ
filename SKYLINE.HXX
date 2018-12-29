//   *********************
//   *** CLASS SKYLINE ***
//   *********************


#ifndef skyline_hxx

#include "column.hxx"
#include <stdio.h>
class Domain ; class FloatMatrix ; class  FloatArray ; class IntArray ;


class Skyline
/*
   This class implements a symmetric matrix stored in a compact (skyline)
   form. A skyline is usually an attribute of the linear system.
 DESCRIPTION :
   Attribute 'size' is number of columns of the skyline. Attribute 'columns'
   is a list of 'size' columns, each of them of any size. Attribute 'isFac-
   torized' is True if the skyline is already in U(T).D.U factorized form,
   else it is False.
 TASKS :
   - storing or returning a coefficient (method 'at').
   - assembling to itself a matrix, typically and elemental stiffness matrix
     (method 'assemble') ; this may include enlarging the profile (method
     'growTo') ;
   - when associated with an array b, solving the system Ax=b (methods 'fac-
     torized','forwardReductionWith',etc) ;
   - when associated with a domain, shaping its own profile (method 'carve-
     Yourself') ;
   - resetting to zero all of its coefficients (method 'reinitialized').
*/
{
   enum { FALSE , TRUE } ;

   protected :
      int           size ;
      Column**      columns ;
      int           isFactorized ;

   public :
      Skyline () ;                         // constructor
      ~Skyline () ;                        // destructor

      FloatMatrix*  AsFloatMatrix () ;
      double&       at (int i,int j)     {return this->giveColumn(j)->at(i);}
      void          assemble (FloatMatrix*,IntArray*) ;
      FloatArray*   backSubstitutionWith (FloatArray*) ;
      void          carveYourselfFor (Domain*) ;
      void          checkSizeTowards (IntArray*) ;
      void          createColumn (int,int) ;
      Skyline*      diagonalScalingWith (FloatArray*) ;
      Skyline*      factorized () ;
      Skyline*      forwardReductionWith (FloatArray*) ;
      Column*       giveColumn (int j)               { return columns[j-1] ;}
      int           giveNumberOfColumns ()           { return size ;}
      void          growTo (int) ;
      void          printYourself () ;
      Skyline*      reinitialized () ;
} ;

#define skyline_hxx
#endif
