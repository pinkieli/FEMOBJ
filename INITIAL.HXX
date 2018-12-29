//   *******************************
//   *** CLASS INITIAL CONDITION ***
//   *******************************
 

#include "load.hxx"
#include "dictionr.hxx"


class InitialCondition : public Load
/*
   This class implements an initial condition. An initial condition is usu-
   ally attribute of one or more degrees of freedom.
*/
{
   private :
      Dictionary*  initialValueDictionary ;

   public :
      InitialCondition (int i,Domain* d) : Load(i,d)
					  { initialValueDictionary = NULL ;}
      ~InitialCondition ()                { delete initialValueDictionary ;}

      double  give (char) ;
      int     hasConditionOn (char u) ;
      void    instanciateYourself ()      { this -> readConditions();}
      void    printYourself () ;
      void    readConditions () ;
} ;

