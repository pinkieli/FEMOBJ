//   ***************************
//   *** CLASS LINEAR SYSTEM ***
//   ***************************
 

#ifndef linsyst_hxx

#include "flotarry.hxx"
#include "skyline.hxx"


class LinearSystem
/*
 DESCRIPTION :
   This class defines the attribute 'leftHandSide', 'rightHandSide' and
   'solutionArray'. Attribute 'currentSize' is used during the equation
   numbering process (see below).
 TASKS :
   - solving itself (implemented by the subclasses of LinearSystem) ;
   - assigning an equation number to any degree of freedom which asks for
     one (method 'giveUpdatedCurrentSize'). The system does so by returning
     its attribute 'currentSize' incremented every times by 1.
*/
{
   protected :
      Skyline*     leftHandSide ;
      FloatArray*  rightHandSide ;
      FloatArray*  solutionArray ;
      int          currentSize ;

   public :
      LinearSystem () ;
      ~LinearSystem () ;

      void          carveYourselfFor (Domain*) ;
      Skyline*      giveLhs ()                    { return leftHandSide ;}
      FloatArray*   giveRhs ()                    { return rightHandSide ;}
      FloatArray*   giveSolutionArray ()          { return solutionArray ;}
      int           giveUpdatedSize ()            { return ++currentSize ;}
      void          solveYourself () ;
      void          solveYourselfAndCheckSolution () ;
      void          updateYourself () ;
      void          updateYourselfExceptLhs () ;
} ;

#define linsyst_hxx
#endif








