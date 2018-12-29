//   file LINSYST.CXX

#include "linsyst.hxx"
#include "flotmtrx.hxx"
#include "flotarry.hxx"
#include "debug.def"
#include "verbose.def"


LinearSystem :: LinearSystem ()
   // Constructor.
{
   leftHandSide  = new Skyline() ;
   rightHandSide = new FloatArray(0) ;
   solutionArray = NULL ;
   currentSize   = 0 ;
}

LinearSystem :: ~LinearSystem ()
   // Destructor.
{
   delete leftHandSide ;
   delete rightHandSide ;
   delete solutionArray ;
}


void  LinearSystem :: carveYourselfFor (Domain* aDomain)
   // Constructs the profile (without computing the coefficients) of the
   // left-hand and the right-hand side of the receiver, associated with the
   // elements and nodes of aDomain.
{
   int n ;

#  ifdef VERBOSE
      printf ("\ncarving the linear system\n") ;
#  endif

   leftHandSide -> carveYourselfFor(aDomain) ;

   n = leftHandSide->giveNumberOfColumns() ;
   delete rightHandSide ;
   rightHandSide = new FloatArray(n) ;
}


void  LinearSystem :: solveYourself ()
   // Solves the system  A.x = b  implemented by the receiver.
{
   FloatArray* y ;         // y is progressively changed into 'solutionArray'

   y = rightHandSide -> GiveCopy() ;
   solutionArray = leftHandSide -> factorized()
				-> forwardReductionWith(y)
				-> diagonalScalingWith (y)
				-> backSubstitutionWith(y) ;
}


void  LinearSystem :: updateYourself ()
   // Updates the receiver at end of time step.
{
   this -> updateYourselfExceptLhs() ;
   leftHandSide  -> reinitialized() ;
}


void  LinearSystem :: updateYourselfExceptLhs ()
   // Updates the receiver at end of time step. Does not reinitialize the
   // left-hand side, which remains unmodified (in a factorized form).
{
   rightHandSide -> reinitialized() ;
   delete solutionArray ;
   solutionArray = NULL ;
}


