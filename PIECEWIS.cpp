//   file PIECEWIS.CXX
 
#include "piecewis.hxx"
#include <math.h>


double  PiecewiseLinFunction :: at (double time)
   // Returns the value of the receiver at time 'time'. 'time' should be
   // one of the dates of the receiver (currently there is no interpola-
   // tion between two points).
{
   const double  precision = 0.000001 ;
   int           i ;

   if (! numberOfPoints)
      this -> getPoints() ;

   for (i=0 ; i<numberOfPoints ; i++)
      if (fabs(dates[i]-time) < precision)
	 return values[i] ;

   return 0. ;
}


void  PiecewiseLinFunction :: getPoints ()
   // Reads in the data file the date and the value of every point of the
   // receiver.
{
   int i ;

   numberOfPoints = this->readInteger("nPoints") ;

   dates = new double [numberOfPoints] ;
   for (i=0 ; i<numberOfPoints ; i++)
      dates[i] = this -> read("t",i+1) ;

   values = new double [numberOfPoints] ;
   for (i=0 ; i<numberOfPoints ; i++)
      values[i] = this -> read("f(t)",i+1) ;
}


