//   file FLOTARRY.CXX

#include "flotarry.hxx"
#include "intarray.hxx"
#include "flotmtrx.hxx"
#include "mathfem.h"
#include <math.h>
#include <string.h>


FloatArray :: FloatArray (int n)
   // Constructor : creates an array of size n (filled with garbage).
{
   size = n ;
   if (size)
      values = allocDouble(size) ;
   else
      values = NULL ;
}


FloatArray*  FloatArray :: add (FloatArray* b)
   // Performs the operation a=a+b, where a stands for the receiver. If the
   // receiver's size is 0, adjusts its size to that of b. Returns the
   // receiver.
{
   register int i ;
   double       *p1,*p2 ;

   if (!b || b->giveSize()==0)
      return this ;

   if (! size) {                              // null-sized array
      size   = b -> size ;
      values = allocDouble(size) ;}

#  ifdef DEBUG
      if (size != b->size) {                  // unmatching sizes
	 printf ("FloatArray dimension mismatch in a[%d]->add(b[%d])\n",
		  size,b->size) ;
	 exit(0) ;}
#  endif

   p1 = values ;
   p2 = b -> values ;
   i  = size ;
   while (i--)
      *p1++ += *p2++ ;
   return this ;
}


void  FloatArray :: assemble (FloatArray* fe, IntArray* loc)
   // Assembles the array fe (typically, the load vector of a finite
   // element) to the receiver, using loc as location array.
{
   int i,ii,n ;

#  ifdef DEBUG
      if ((n=fe->giveSize()) != loc->giveSize()) {
	 printf ("dimensions of 'fe' and 'loc' mismatch \n") ;
	 exit(0) ; }
      this -> checkSizeTowards(loc) ;
#  endif

   n = fe->giveSize() ;
   for (i=1 ; i<=n ; i++) {
      ii = loc->at(i) ;
      if (ii)                                  // if non 0 coefficient,
	 this->at(ii) += fe->at(i) ; }         // then assemble
}


#ifdef DEBUG
double&  FloatArray :: at (int i)
   // Returns the i-th coefficient of the receiver. Slow but safe.
{
   this -> checkBounds(i) ;
   return values[i-1] ;
}
#endif


#ifdef DEBUG
void  FloatArray :: checkBounds(int i)
   // Checks that the receiver's size is not smaller than 'i'.
{
   if (i<=0) {
      printf ("array error on index : %d <= 0 \n",i) ;
      exit(0) ; }

   if (i>size) {
      printf ("array error on index : %d > %d \n",i,size) ;
      exit(0) ; }
}
#endif


void  FloatArray :: checkSizeTowards (IntArray* loc)
   // Expands the receiver if loc points to coefficients beyond the size of
   // the receiver.
{
   int i,n,high ;

   high = 0 ;
   n    = loc -> giveSize() ;
   for (i=1 ; i<=n ; i++)
      high = max (high,(loc->at(i))) ;
   if (high > size)                             // receiver must be expanded
      this -> growTo(high) ;
}


int  FloatArray :: containsOnlyZeroes ()
   // Returns True if all coefficients of the receiver are 0, else returns
   // False.
{
   register int i ;
   double       *p ;

   p = values ;
   i = size ;
   while (i--)
      if (*p++ != 0.)
	 return FALSE ;

   return TRUE ;
}


void  FloatArray :: growTo (int n)
   // Expands the receiver up to size n (n is assumed larger than 'size').
   // Initializes all new coefficients to zero.
{
   register int i ;
   double       *newValues,*p1,*p2 ;

#  ifdef DEBUG
      if (!n || n<=size) {
	 printf ("error in FloatArray(%d)::growTo(%d) \n",size,n) ;
	 exit(0) ;}
#  endif

   newValues = allocDouble(n) ;

   p1 = values ;
   p2 = newValues ;
   i  = size ;
   while (i--)
      *p2++ = *p1++ ;

   if (values)
      freeDouble (values) ;
   values = newValues ;
   size   = n ;
}


FloatArray*  FloatArray :: negated ()
   // Switches the sign of every coefficient of the receiver. Returns the
   // receiver.
{
   register int i ;
   double       x ;
   double*      p ;

   i = size ;
   p = values ;
   while (i--) {
      x    = - *p ;
      *p++ = x ;}
   return this ;
}


void  FloatArray :: printYourself ()
   // Prints the receiver on screen.
{
   printf ("FloatArray of size : %d \n",size) ;
   for (int i=1 ; i<=size ; ++i)
      printf ("%.12f  ",this->at(i)) ;
   printf ("\n") ;
}


FloatArray*  FloatArray :: reinitialized ()
   // Returns the receiver with all coefficients set to 0.
{
   if (values)
      freeDouble (values) ;
   values = allocDouble(size) ;
   return this ;
}


FloatArray*  FloatArray :: rotatedWith (FloatMatrix* r, char mode)
   // Returns the receiver 'a' rotated according the change-of-base matrix r.
   // If mode = 't', the method performs the operation  a = r(transp) * a .
   // If mode = 'n', the method performs the operation  a = r * a .
{
   register int i ;
   double       *p1,*p2 ;
   FloatMatrix  *rot ;
   FloatArray   *rta ;

   if (mode == 't')
      rot = r -> GiveTransposition() ;
   else
      rot = r ;

   rta = rot -> Times(this) ;

   p1 = values ;
   p2 = rta -> values ;
   i  = size ;
   while (i--)
      *p1++ = *p2++ ;

   if (mode == 't')
      delete rot ;
   delete rta ;
   return this ;
}


FloatArray*  FloatArray :: times (double factor)
   // Multiplies every coefficient of the receiver by factor. Answers the
   // modified receiver.
{
   register int i ;
   double*      p ;

   p = values ;
   i = size ;
   while (i--)
      *(p++) *= factor ;
   return this ;
}


FloatArray*  FloatArray :: Times (double factor)
   // Returns a new array, whose components are those of the receicer, times
   // factor.
{
   register int i ;
   double       *p1,*p2 ;
   FloatArray*  answer ;

   answer = new FloatArray(size) ;
   p1     = values ;
   p2     = answer -> values ;
   i      = size ;
   while (i--)
      *p2++ = factor * (*p1++) ;
   return answer ;
}


double  dotProduct (double* P1, double* P2, register int i)
   // A non-member function. Returns the dot product of the first 'i' coef-
   // ficienst of the two arrays P1 and P2. This method applies to many
   // situations, eg row*column products with matrices.
{
   double answer ;

   answer = 0. ;
   while (i--)
      answer += *P1++ * *P2++ ;
   return answer ;
}

