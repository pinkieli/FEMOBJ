//   file NODLOAD.CXX
 
#include "nodload.hxx"
#include "loadtime.hxx"
#include "flotarry.hxx"
#include "timestep.hxx"


FloatArray*  NodalLoad :: ComputeValueAt (TimeStep* stepN)
   // Returns an array, the force induced at stepN by the receiver.
{
   FloatArray* force ;
   double      factor ;

   factor = this -> giveLoadTimeFunction() -> at(stepN->giveTime()) ;
   if (factor == 0.)
      return NULL ;
   else {
      force  = this -> giveComponentArray() ;
      return  force -> Times(factor) ;}
}







