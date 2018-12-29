//   file NODE.CXX

#include "node.hxx"
#include "dof.hxx"
#include "nodload.hxx"
#include "timestep.hxx"
#include "linsyst.hxx"
#include "flotarry.hxx"
#include "intarray.hxx"
#include "debug.def"
#include "verbose.def"
#include <stdlib.h>


Node :: Node (int n, Domain* aDomain)
      : FEMComponent (n,aDomain)
   // Constructor. Creates a node with number n, belonging to aDomain.
{
   numberOfDofs  = 0 ;
   coordinates   = NULL ;
   dofArray      = NULL ;
   loadArray     = NULL ;
   locationArray = NULL ;
}


Node :: ~Node()
   // Destructor.
{
   int i = numberOfDofs ;

   delete coordinates ;
   if (numberOfDofs) {
      while (i--)
	 delete dofArray[i] ;
      //delete [numberOfDofs] dofArray ; // revised by pinkie.ljz
      delete [] dofArray ;
      }
   delete loadArray ;
   delete locationArray ;
}


void  Node :: assembleYourLoadsAt (TimeStep* stepN)
   // Forms at stepN the vector of the concentrated loads of the receiver ;
   // then, if it exists, assembles it to the system's right-hand side.
{
   LinearSystem *system ;
   FloatArray   *loadVector,*rhs ;
   IntArray     *loc ;

#  ifdef VERBOSE
      printf ("assembling loads of node %d\n",number) ;
#  endif

   loadVector = this -> ComputeLoadVectorAt(stepN) ;
   if (loadVector) {
      system = (LinearSystem*) (domain->giveLinearSystem()) ;
      rhs    = system -> giveRhs() ;
      loc    = this -> giveLocationArray() ;
      rhs -> assemble(loadVector,loc) ;
      delete loadVector ;}
}


FloatArray*  Node :: ComputeLoadVectorAt (TimeStep* stepN)
   // Computes the vector of the nodal loads of the receiver.
{
   int        i,n,nLoads ;
   NodalLoad  *loadN ;
   FloatArray *answer,*contribution ;

   if (this -> giveLoadArray() -> isEmpty())
      return NULL ;

   else {
      answer = new FloatArray(0) ;
      nLoads = loadArray->giveSize() ;          // the node may be subjected
      for (i=1 ; i<=nLoads ; i++) {             // to more than one load
	 n            = loadArray -> at(i) ;
	 loadN        = (NodalLoad*) domain->giveLoad(n) ;
	 contribution = loadN -> ComputeValueAt(stepN) ;      // can be NULL
	 answer -> add(contribution) ;
	 delete contribution ;}
      if (answer->giveSize())
	 return answer ;
      else {
	 delete answer ;
	 return NULL ;}}
}


void  Node :: getCoordinates ()
   // Get from the data file all of the coordinates of the receiver.
{
   int  i,numberOfCoordinates ;

   numberOfCoordinates = this->readInteger("coord") ;
   coordinates         = new FloatArray(numberOfCoordinates) ;
   for (i=1 ; i<=numberOfCoordinates ; i++)
      coordinates->at(i) = this->read("coord",i+1) ;
}


double  Node :: giveCoordinate (int i)
   // Returns the i-th coordinate of the receiver.
{
   if (! coordinates)
      this -> getCoordinates() ;
   return coordinates->at(i) ;
}


Dof*  Node :: giveDof (int i)
   // Returns the i-th degree of freedom of the receiver. Creates the array
   // containing the dofs of the receiver, if such array does not exist yet.
{
   int j ;

   if (! dofArray) {
      dofArray = new Dof* [this->giveNumberOfDofs()] ;
      for (j=0 ; j<numberOfDofs ; j++)
	 dofArray[j] = new Dof(j+1,this) ;}

   return dofArray[i-1] ;
}


IntArray*  Node :: giveLoadArray ()
   // Returns the list containing the number of every nodal loads that act on
   // the receiver. If this list does not exist yet, constructs it. This list
   // is not to be confused with the load vector.
{
   int i,nLoads ;

   if (! loadArray) {                          // the list does not exist yet
      nLoads    = this->readIfHas("loads") ;
      loadArray = new IntArray(nLoads) ;
      for (i=1 ; i<=nLoads ; i++)
	 loadArray->at(i) = this->readInteger("loads",i+1) ;}

   return loadArray ;
}


IntArray*  Node :: giveLocationArray ()
   // Returns the location array of the receiver. Creates this array if it
   // does not exist yet. The location array contains the equation number of
   // every  degree of freedom of the receiver.
{
   int i ;

   if (! locationArray) {
      locationArray = new IntArray(this->giveNumberOfDofs()) ;
      for (i=1 ; i<=numberOfDofs ; i++)
	 locationArray->at(i) = this->giveDof(i)->giveEquationNumber() ;}

   return locationArray ;
}


int  Node :: giveNumberOfDofs ()
   // Returns the number of degrees of freedom of the receiver.
{
   if (! numberOfDofs)
      numberOfDofs = this->readInteger("nDofs") ;

   return numberOfDofs ;
}


void  Node :: instanciateYourself ()
   // Gets from the data file all the data of the receiver.
{
   int i ;

#  ifdef VERBOSE
      printf ("instanciating node %d\n",number) ;
#  endif

   this -> getCoordinates() ;
   this -> giveLocationArray () ;
   this -> giveLoadArray() ;

   numberOfDofs = this -> readInteger("nDofs") ;
   for (i=1 ; i<=numberOfDofs ; i++)
      this -> giveDof(i) -> hasBc() ;
}


void  Node :: printOutputAt (TimeStep* stepN)
{
   int  i ;

#  ifdef VERBOSE
      printf ("node %d printing output\n",number) ;
#  endif

   for (i=1 ; i<=numberOfDofs ; i++)
      this -> giveDof(i) -> printOutputAt(stepN) ;
}


void  Node :: printYourself ()
   // Prints the receiver on screen.
{
   int    i ;
   double x,y ;

   x = this->giveCoordinate(1) ;
   y = this->giveCoordinate(2) ;
   printf ("Node %d    coord : x %f  y %f\n",number,x,y) ;
   for (i=0 ; i<numberOfDofs ; i++) {
      if (dofArray[i])
	 dofArray[i] -> printYourself() ;
      else
	 printf ("dof %d is nil \n",i+1) ;}
   if (locationArray)
      locationArray->printYourself() ;
   else
      printf ("locationArray = nil \n") ;
   if (loadArray)
      loadArray->printYourself() ;
   else
      printf ("loadArray = nil \n") ;
   printf ("\n") ;
}


void  Node :: updateYourself ()
   // Updates the receiver at end of step.
{
   int i ;

#  ifdef VERBOSE
      printf ("updating node %d\n",number) ;
#  endif

   for (i=1 ; i<=numberOfDofs ; i++)
      this -> giveDof(i) -> updateYourself() ;
}
