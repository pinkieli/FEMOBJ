//   file DEADWGHT.CXX
 
#include "deadwght.hxx"
#include "element.hxx"
#include "timestep.hxx"
#include "loadtime.hxx"
#include "material.hxx"
#include "flotarry.hxx"


FloatArray*  DeadWeight :: ComputeForceOn (Element* elem, TimeStep* stepN)
   // Returns the force (and not just the acceleration) induced at stepN by
   // the receiver on the element 'elem'.
{
   double density,factor ;

   factor = this -> giveLoadTimeFunction() -> at(stepN->giveTime()) ;
   if (factor == 0.)
      return NULL ;
   else {
      density = elem -> giveMaterial() -> give('d') ;
      return  this -> giveComponentArray() -> Times(density*factor) ;}
}


