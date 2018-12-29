//   ***********************
//   *** CLASS BODY LOAD ***
//   ***********************


#ifndef bodyload_hxx

#include "load.hxx"
class Element ; class TimeStep ;


class BodyLoad : public Load
/*
   This abstract class is the superclass of all loads (e.g., the dead weight)
   that act on the whole volume of finite elements. A body load is usually
   attribute of many elements.
*/

{
   public :
      BodyLoad (int i,Domain* d) : Load (i,d) {}        // constructor

      virtual FloatArray*  ComputeForceOn (Element*,TimeStep*) = 0 ; // revised by pinkie.ljz; null is replaced by 0.
} ;

#define bodyload_hxx
#endif








