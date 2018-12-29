//   ******************
//   *** CLASS NODE ***
//   ******************


#include "femcmpnn.hxx"
#include "domain.hxx"
#include <stdio.h>
class Dof ; class NodalLoad ; class TimeStep ;
class FloatArray ; class IntArray ;


class Node : public FEMComponent
/*
   This class implements a node in a finite element mesh. A node is an attri-
   bute of a domain. It is usually also attribute of a few elements.
 DESCRIPTION
   The node possesses 'numberOfDofs' degrees of freedom, stored in 'dofArray'.
   In 'loadArray' it stores the number of every nodal load it is subjected to
   (typically, concentrated forces and moments).
   In 'locationArray' the node stores the equation number of each of its dofs.
   This location array is used by the node for assembling its load vector to
   the right-hand side of the linear system ; it is also used by elements for
   calculating their own location arrays.
 TASKS
   - managing its position in space. In geometrically linear analysis, this
     only amounts to managing its coordinates ;
   - managing its degrees of freedom (method 'giveDof') ;
   - calculating its nodal load vector and having it assembled to the linear
     system ;
   - printing and updating at end of step ;
   - managing its swapping to and from disk.
*/
{
   private:
      FloatArray*  coordinates ;
      int          numberOfDofs ;
      Dof**        dofArray ;
      IntArray*    loadArray ;
      IntArray*    locationArray ;

   public:
      Node (int,Domain*) ;                        // constructor
      ~Node () ;                                  // destructor

      // coordinates
      double       giveCoordinate (int) ;
      void         getCoordinates () ;

      // degrees of freedom
      Dof*         giveDof (int) ;
      int          giveNumberOfDofs () ;

      // nodal load vector
      void         assembleYourLoadsAt (TimeStep*) ;
      FloatArray*  ComputeLoadVectorAt (TimeStep*) ;
      IntArray*    giveLoadArray () ;
      IntArray*    giveLocationArray () ;

      // time step termination
      void         printOutputAt (TimeStep*) ;
      void         updateYourself () ;

      // miscellaneous
      char*        giveClassName (char* s)      { return strcpy(s,"Node") ;}
      void         instanciateYourself () ;
      void         printYourself () ;
} ;









