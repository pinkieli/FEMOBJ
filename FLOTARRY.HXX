//   *************************
//   *** CLASS FLOAT ARRAY ***
//   *************************
 

#ifndef flotarry_hxx

#include "freestor.h"
#include "debug.def"
#include <stdlib.h>
class IntArray ; class FloatMatrix ;


class FloatArray
/*
   This class implements an array of floating-point numbers.
 DESCRIPTION :
   A FloatArray stores its coefficients in an array 'values' of size 'size'.
 TASKS :
   - storing and returning a coefficient (method 'at') ;
   - expanding its size in order to store additional coefficients (method
     'growTo') ;
   - performing basic operations : summation, product, rotation, etc ;
   - assembling to itself another array, typically an elemental or nodal
     load vector (method 'assemble') ;
   - reading/writing its description on a given file.
 REMARKS :
   - for the sake of efficiency, the array 'values' is allocated using the
     C 'calloc' function rather than the 'new' operator.
   - method 'givePointer' is an encapsulation crime. It is used only for
     speeding up method 'dot' of class RowColumn.
*/
{
   enum { FALSE , TRUE } ;

   protected :
      int      size ;
      double*  values ;

   public :
      FloatArray (int) ;                                       // constructor
      virtual ~FloatArray () {if(values) freeDouble(values);}  // destructor

#ifdef DEBUG
      double&      at (int i) ;
#else
      double&      at (int i)                { return values[i-1] ;}
#endif
      FloatArray*  add (FloatArray*) ;
      void         assemble (FloatArray*,IntArray*) ;
      void         checkBounds (int) ;
      void         checkSizeTowards (IntArray*) ;
      int          containsOnlyZeroes () ;
      FloatArray*  GiveCopy ()               { return this->Times(1.) ;}
      double*      givePointer ()            { return values ;} // see above
      int          giveSize ()               { return size ;}
      void         growTo (int) ;
      int          isNotEmpty ()             { return (size != 0) ;}
      FloatArray*  negated () ;
      void         printYourself () ;
      FloatArray*  reinitialized() ;
      FloatArray*  rotatedWith (FloatMatrix*,char) ;
      FloatArray*  times (double) ;
      FloatArray*  Times (double) ;

      friend double  dotProduct (double*,double*,int) ;
} ;

#define flotarry_hxx
#endif








