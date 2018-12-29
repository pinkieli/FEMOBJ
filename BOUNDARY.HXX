//   ********************************
//   *** CLASS BOUNDARY CONDITION ***
//   ********************************
 

#include "load.hxx"
#include "dictionr.hxx"
class TimeStep ; class Dictionary ;


class BoundaryCondition : public Load
/*
   This class implements a kinematic boundary condition. A b.c. is usually
   attribute of one or more degrees of freedom.
 DESCRIPTION
   The inherited attribute 'componentArray is not used. It is replaced with
   the more adequate dictionary 'prescribedValueDictionary', which entries are
   referenced by names rather than by indices.
 TASKS
   returning a component, i.e., the prescribed value of a kinematic unknown
   (displacement, velocity, temperature, etc).
   REMARK
   Like the other Loads, a b.c. possesses a load-time function, which is ty-
   pically a ConstantFunction.
*/
{
   private :
      Dictionary*  prescribedValueDictionary ;

   public :
      BoundaryCondition (int i,Domain* d) : Load(i,d)
				       { prescribedValueDictionary = NULL ;}
      ~BoundaryCondition ()            { delete prescribedValueDictionary ;}

      double  give (char,TimeStep*) ;
      void    instanciateYourself ()   { this->readConditions() ;}
      void    readConditions () ;
} ;

