//   file DOF.CXX
 
#include "dof.hxx"
#include "node.hxx"
#include "domain.hxx"
#include "timestep.hxx"
#include "timinteg.hxx"
#include "boundary.hxx"
#include "initial.hxx"
#include "linsyst.hxx"
#include "flotarry.hxx"
#include "dictionr.hxx"
#include "freader.hxx"
#include "string.hxx"
#include "debug.def"
#include <stdlib.h>
#include <ctype.h>


Dof :: Dof (int i, Node* aNode)
   // Constructor. Creates a new d.o.f., with number i, belonging
   // to aNode.
{
   number         = i ;
   node           = aNode ;
   equationNumber = -1 ;                         // means "uninitialized"
   bc             = -1 ;
   ic             = -1 ;
   unknowns       = new Dictionary() ;           // unknown size ?
   pastUnknowns   = NULL ;
}


double  Dof :: computeNewmarkUnknown (char u, TimeStep* stepN)
   // Returns the value of the unknown 'u' (eg, the displacement 'd') of the
   // receiver, at stepN, using Newmark's formula.
{
   double                 d,v,a,dPred,vPred,dt,beta,gamma ;
   TimeIntegrationScheme* scheme ;

   if (u == 'd') {                                  // displacement
      scheme = node -> giveDomain() -> giveTimeIntegrationScheme() ;
      dPred  = this -> giveUnknown('D',stepN) ;
      a      = this -> giveUnknown('a',stepN) ;
      beta   = scheme -> giveBeta() ;
      dt     = stepN -> giveTimeIncrement() ;
      return  dPred + beta*dt*dt*a ;}

   else if (u == 'v') {                             // velocity
      scheme = node -> giveDomain() -> giveTimeIntegrationScheme() ;
      vPred  = this -> giveUnknown('V',stepN) ;
      a      = this -> giveUnknown('a',stepN) ;
      gamma  = scheme -> giveGamma() ;
      dt     = stepN -> giveTimeIncrement() ;
      return  vPred + gamma*dt*a ;}

   else if (u == 'a')                               // acceleration
      return  this -> getSolution() ;

   else if (u == 'D') {                             // predicted displacement
      scheme = node -> giveDomain() -> giveTimeIntegrationScheme() ;
      d      = pastUnknowns -> at('d') ;
      v      = pastUnknowns -> at('v') ;
      a      = pastUnknowns -> at('a') ;
      beta   = scheme -> giveBeta() ;
      dt     = stepN -> giveTimeIncrement() ;
      return  d + dt*v + (0.5-beta)*dt*dt*a ;}

   else if (u == 'V') {                              // predicted velocity
      scheme = node -> giveDomain() -> giveTimeIntegrationScheme() ;
      v      = pastUnknowns -> at('v') ;
      a      = pastUnknowns -> at('a') ;
      gamma  = scheme -> giveGamma() ;
      dt     = stepN -> giveTimeIncrement() ;
      return  v + (1.-gamma)*dt*a ;}

   else {
      printf ("Error : no such unknown with Newmark : %c\n",u) ;
      exit(0) ;}
   return 0.;
}


double  Dof :: computeStaticUnknown (char u, TimeStep* stepN)
   // Returns the value of the unknown 'u' (necessarily, the displacement
   // 'd') of the receiver, at stepN.
{
   if (u == 'd')                                   // displacement
      return  this -> getSolution() ;
   else {
      printf ("Error : no such unknown with Static : %c\n",u) ;
      exit(0) ;}
   return 0.;
}


double  Dof :: computeUnknown (char u, TimeStep* stepN)
   // Computes the value of the unknown 'u' (e.g., the displacement 'd') of
   // the receiver, at stepN.
{
   TimeIntegrationScheme* scheme ;

   scheme = node -> giveDomain() -> giveTimeIntegrationScheme() ;
   if (scheme -> isStatic())
      return  this -> computeStaticUnknown(u,stepN) ;
   else if (scheme->isNewmark())
      return  this -> computeNewmarkUnknown(u,stepN) ;
   else {
      printf ("Error : unknown time integration scheme : %c\n",scheme) ;
      exit(0) ;}
   return 0.;
}


double  Dof :: getSolution ()
   // Gets from the linear system the solution associated to the receiver.
{
   LinearSystem* System ;
   double        answer ;

   System = (LinearSystem*) (node->giveDomain()->giveLinearSystem()) ;
   answer = System -> giveSolutionArray() -> at(equationNumber) ;
   return answer ;
}


BoundaryCondition*  Dof :: giveBc ()
   // Returns the boundary condition the receiver is subjected to.
{
#  ifdef DEBUG
      if (bc == -1) {
	 printf ("error : Dof does not know yet if has bc or not \n") ;
	 exit(0) ;}
#  endif

   return  (BoundaryCondition*) (node -> giveDomain() -> giveLoad(bc)) ;
}


int  Dof :: giveEquationNumber ()
   // Returns the number of the equation in the linear system that corres-
   // ponds to the receiver. The equation number is 0 if the receiver is
   // subjected to a boundary condition, else it is n+1, where n is the
   // equation number of the most recently numbered degree of freedom.
{
   LinearSystem* system ;

   if (equationNumber == -1) {                      // not yet computed
      if (this -> hasBc())
	 equationNumber = 0 ;
      else {
	 system = (LinearSystem*) (node->giveDomain()->giveLinearSystem()) ;
	 equationNumber = system -> giveUpdatedSize() ;}}

   return equationNumber ;
}


InitialCondition*  Dof :: giveIc ()
   // Returns the initial condition on the receiver. Not used.
{
#  ifdef DEBUG
      if (ic == -1) {
	 printf ("error : Dof does not know yet if has InitCond or not \n") ;
	 exit(0) ;}
#  endif

   return  (InitialCondition*) (node -> giveDomain() -> giveLoad(ic)) ;
}


double  Dof :: givePastUnknown (char u, TimeStep* stepN)
   // Returns the value of the unknown 'u' (e.g., the displacement) of the
   // receiver, at a previous time step stepN. If stepN is step 0, the value
   // is given by the initial condition, if any.
{
   double value ;

/*
   if (stepN->giveNumber() == 1) {
      if (pastUnknowns -> includes(u))
	 value = pastUnknowns -> at(u) ;
      else {
	 if (this -> hasIc())
	    value = this -> giveIc() -> give(u) ;
	 else
	    value = 0. ;
	 pastUnknowns -> add(u,value) ;}}

   else
*/
      value = pastUnknowns -> at(u) ;

   return value ;
}


double  Dof :: giveUnknown (char u, TimeStep* stepN)
   // The key method of class Dof. Returns the value of the unknown 'u'
   // (e.g., the displacement) of the receiver, at stepN. This value may,
   // or may not, be already available. It may depend on a boundary (if it
   // is not a predicted unknown) or initial condition. stepN is not the
   // current time step n, it is assumed to be the previous one (n-1).
{
   double value ;

   if (stepN -> isTheCurrentTimeStep()) {
      if (unknowns -> includes(u))                       // already known
	 value = unknowns -> at(u) ;
      else {
	 if (stepN->giveNumber()==0) {                   // step 0
	    if (this->hasIcOn(u))                        //   init. cond.
	       value = this -> giveIc() -> give(u) ;
	    else                                         //   no init. cond.
	       value = 0. ;}
	 else if (this->hasBc() && islower(u))           // bound . cond.
	    value = this -> giveBc() -> give(u,stepN) ;
	 else                                            // compute it !
	    value = this -> computeUnknown(u,stepN) ;
	 unknowns -> add(u,value) ;}}

   else
      value = this -> givePastUnknown(u,stepN) ;         // the previous step

   return value ;
}



int  Dof :: hasBc ()
   // Returns True if the receiver is subjected to a boundary condition, else
   // returns False. If necessary, reads the answer in the data file.
{
   char bcOnDofN[32] ;

   if (bc == -1) {
      concatenate("bcOnDof",number,bcOnDofN) ;
      bc = node -> readIfHas(bcOnDofN) ;}

   return bc ;
}


int  Dof :: hasIc ()
   // Returns True if the receiver is subjected to an initial condition,
   // else returns False.
{
   char icOnDofN[32] ;

   if (ic == -1) {
      concatenate("icOnDof",number,icOnDofN) ;
      ic = node -> readIfHas(icOnDofN) ;}

   return ic ;
}


int  Dof :: hasIcOn (char u)
   // Returns True if the unknown 'u' (e.g., the displacement 'd') of the
   // receiver is subjected to an initial condition, else returns False.
{
   if (this->hasIc())
      return this->giveIc()->hasConditionOn(u) ;
   else
      return FALSE ;
}


void  Dof :: print (char u, TimeStep* stepN)
   // Prints in the data file the unknown 'u' (for example, the displacement
   // 'd') of the receiver, at stepN.
{
   FILE*  File ;
   int    nod ;
   double x ;

   File = node -> giveDomain() -> giveOutputStream() ;
   nod  = node -> giveNumber() ;
   x    = this -> giveUnknown(u,stepN) ;
   if (number == 1)
      fprintf (File,"node%4d  ",nod) ;
   else
      fprintf (File,"          ") ;
   fprintf (File,"dof %d   %c % .8e\n",number,u,x) ;
}


void  Dof :: print (char u1,char u2,char u3,TimeStep* stepN)
   // Prints in the data file the unknowns 'u1' (for example, the displacement
   // 'd'), 'u2' and 'u3' of the receiver, at stepN.
{
   FILE*  File ;
   int    nod ;
   double x1,x2,x3 ;

   File = node -> giveDomain() -> giveOutputStream() ;
   nod  = node -> giveNumber() ;
   x1   = this -> giveUnknown(u1,stepN) ;
   x2   = this -> giveUnknown(u2,stepN) ;
   x3   = this -> giveUnknown(u3,stepN) ;
   if (number == 1)
      fprintf (File,"node%4d  ",nod) ;
   else
      fprintf (File,"          ") ;
   fprintf (File,"dof %d   %c % .8e   %c % .8e   %c % .8e\n",
      number,u1,x1,u2,x2,u3,x3) ;
}


void  Dof :: printOutputAt (TimeStep* stepN)
   // Prints in the data file the unknowns of the receiver at time step
   // stepN. Switches to the correct formula.
{
   TimeIntegrationScheme* scheme ;

   scheme = node -> giveDomain() -> giveTimeIntegrationScheme() ;
   if (scheme -> isStatic())
      this -> printStaticOutputAt(stepN) ;
   else if (scheme->isNewmark())
      this -> printNewmarkOutputAt(stepN) ;
   else {
      printf ("Error : unknown time integration scheme : %c\n",scheme) ;
      exit(0) ;}
}


void  Dof :: printYourself ()
   // Prints the receiver on screen.
{
   printf ("dof %d  of node %d :\n",number,node->giveNumber()) ;
   printf ("equation %d    bc %d \n",equationNumber,bc) ;

   printf ("unknowns : ") ;
   unknowns->printYourself() ;
   if (pastUnknowns) {
      printf ("pastUnknowns : ") ;
      pastUnknowns -> printYourself() ;}
}


void  Dof :: updateYourself ()
   // Updates the receiver at end of step.
{
   delete pastUnknowns ;
   pastUnknowns = unknowns ;
   unknowns     = new Dictionary() ;
}




