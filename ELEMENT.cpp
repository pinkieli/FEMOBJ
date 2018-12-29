//   file ELEMENT.CXX

#include "element.hxx"
#include "planstrn.hxx"
#include "truss2d.hxx"
#include "domain.hxx"
#include "timestep.hxx"
#include "timinteg.hxx"
#include "node.hxx"
#include "dof.hxx"
#include "material.hxx"
#include "bodyload.hxx"
#include "gausspnt.hxx"
#include "intarray.hxx"
#include "flotarry.hxx"
#include "flotmtrx.hxx"
#include "linsyst.hxx"
#include "skyline.hxx"
#include "debug.def"
#include "verbose.def"
#include <stdlib.h>
#include <stdio.h>


Element :: Element (int n, Domain* aDomain)
	 : FEMComponent (n, aDomain)
   // Constructor. Creates an element with number n, belonging to aDomain.
{
   material           = 0    ;
   numberOfNodes      = 0    ;
   nodeArray          = NULL ;
   locationArray      = NULL ;
   constitutiveMatrix = NULL ;
   massMatrix         = NULL ;
   stiffnessMatrix    = NULL ;
   bodyLoadArray      = NULL ;
   gaussPointArray    = NULL ;
}


Element :: ~Element ()
   // Destructor.
{
   int i ;

   delete nodeArray ;
   delete locationArray ;
   delete massMatrix ;
   delete stiffnessMatrix ;
   delete constitutiveMatrix ;
   if (gaussPointArray) {
      for (i=0 ; i<numberOfGaussPoints ; i++)
	 delete gaussPointArray[i] ;
      //delete [numberOfGaussPoints] gaussPointArray ; // revised by pinkie.ljz
      delete [] gaussPointArray ;
      }
   delete bodyLoadArray ;
}


void  Element :: assembleLhsAt (TimeStep* stepN)
   // Assembles the left-hand side (stiffness matrix) of the receiver to
   // the linear system' left-hand side, at stepN.
{
   FloatMatrix* elemLhs ;
   Skyline*     systLhs ;
   IntArray*    locArray ;

   elemLhs  = this -> ComputeLhsAt(stepN) ;
   systLhs  = domain -> giveLinearSystem() -> giveLhs() ;
   locArray = this -> giveLocationArray() ;
   systLhs -> assemble(elemLhs,locArray) ;

   delete elemLhs ;
}


void  Element :: assembleRhsAt (TimeStep* stepN)
   // Assembles the right-hand side (load vector) of the receiver to
   // the linear system' right-hand side, at stepN.
{
   FloatArray* elemRhs ;
   FloatArray* systRhs ;
   IntArray*   locArray ;

   elemRhs = this -> ComputeRhsAt(stepN) ;
   if (elemRhs) {
      systRhs  = domain -> giveLinearSystem() -> giveRhs() ;
      locArray = this -> giveLocationArray() ;
      systRhs -> assemble(elemRhs,locArray) ;
      delete elemRhs ;}
}


void  Element :: assembleYourselfAt (TimeStep* stepN)
   // Assembles the contributions of the receiver to the linear system, at
   // time step stepN. This may, or may not, require assembling the receiver's
   // left-hand side.
{
#  ifdef VERBOSE
      printf ("assembling element %d\n",number) ;
#  endif

   if (stepN -> requiresNewLhs())
      this -> assembleLhsAt(stepN) ;
   this -> assembleRhsAt(stepN) ;
}


FloatArray*  Element :: ComputeBcLoadVectorAt (TimeStep* stepN)
   // Computes the load vector due to the boundary conditions acting on the
   // receiver's nodes, at stepN. Returns NULL if this array contains only
   // zeroes.
{
   FloatArray *d,*answer ;

   d = this -> ComputeVectorOfPrescribed('d',stepN) ;
   if (d -> containsOnlyZeroes())
      answer = NULL ;
   else
      answer = this -> giveStiffnessMatrix() -> Times(d) -> negated() ;

   delete d ;
   return answer ;
}


FloatArray*  Element :: ComputeBodyLoadVectorAt (TimeStep* stepN)
   // Computes numerically the load vector of the receiver due to the body
   // loads, at stepN.
{
   int         i ;
   double      dV ;
   GaussPoint* gp ;
   FloatArray  *answer,*f,*ntf ;
   FloatMatrix *n,*nt ;

   if (this -> giveBodyLoadArray() -> isEmpty())         // no loads
      return NULL ;

   else {
      f = this -> ComputeResultingBodyForceAt(stepN) ;
      if (! f)                                           // nil resultant
	 return NULL ;
      else {
	 answer = new FloatArray(0) ;
	 for (i=0 ; i<numberOfGaussPoints ; i++) {
	    gp  = gaussPointArray[i] ;
	    n   = this -> ComputeNmatrixAt(gp) ;
	    dV  = this -> computeVolumeAround(gp) ;
	    nt  = n    -> GiveTransposition() ;
	    ntf = nt   -> Times(f) -> times(dV) ;
	    answer -> add(ntf) ;
	    delete n ;
	    delete nt ;
	    delete ntf ;}
	 delete f ;
	 return answer ;}}
}


FloatMatrix*  Element :: ComputeConsistentMassMatrix ()
   // Computes numerically the consistent (full) mass matrix of the receiver.
{
   int         i ;
   double      density,dV ;
   FloatMatrix *n,*answer ;
   GaussPoint  *gp ;

   answer  = new FloatMatrix() ;
   density = this -> giveMaterial() -> give('d') ;
   for (i=0 ; i<numberOfGaussPoints ; i++) {
      gp      = gaussPointArray[i] ;
      n       = this -> ComputeNmatrixAt(gp) ;
      dV      = this -> computeVolumeAround(gp) ;
      answer -> plusProduct(n,n,density*dV) ;
      delete n ;}

   return  answer->symmetrized() ;
}


FloatMatrix*  Element :: ComputeLhsAt (TimeStep* stepN)
   // Computes the contribution of the receiver to the left-hand side of the
   // linear system.
{
   TimeIntegrationScheme* scheme ;

   scheme = domain -> giveTimeIntegrationScheme() ;
   if (scheme -> isStatic())
      return  this -> ComputeStaticLhsAt (stepN) ;
   else if (scheme -> isNewmark())
      return  this -> ComputeNewmarkLhsAt(stepN) ;
   else {
      printf ("Error : unknown time integration scheme : %c\n",scheme) ;
      exit(0) ;}
   return new FloatMatrix();
}


FloatArray*  Element :: ComputeLoadVectorAt (TimeStep* stepN)
   // Computes the load vector of the receiver, at stepN.
{
   FloatArray* loadVector ;
   FloatArray* bodyLoadVector = NULL ;
   FloatArray* bcLoadVector   = NULL ;

   loadVector = new FloatArray(0) ;

   bodyLoadVector = this -> ComputeBodyLoadVectorAt(stepN) ;
   if (bodyLoadVector) {
      loadVector -> add(bodyLoadVector) ;
      delete bodyLoadVector ;}

   bcLoadVector = this -> ComputeBcLoadVectorAt(stepN) ;
   if (bcLoadVector) {
      loadVector -> add(bcLoadVector) ;
      delete bcLoadVector ;}

   if (loadVector -> isNotEmpty())
      return loadVector ;
   else {
      delete loadVector ;
      return NULL ;}
}


FloatMatrix*  Element :: computeMassMatrix ()
   // Returns the lumped mass matrix of the receiver.
{
   FloatMatrix* consistentMatrix ;

   consistentMatrix = this -> ComputeConsistentMassMatrix() ;
   massMatrix       = consistentMatrix -> Lumped() ;
   delete consistentMatrix ;

   return massMatrix ;
}


FloatMatrix*  Element :: ComputeNewmarkLhsAt (TimeStep* stepN)
   // Computes the contribution of the receiver to the left-hand side of the
   // linear system, using Newmark's formula.
{
   FloatMatrix *m,*k,*lhs ;
   double      beta,dt ;

   if (stepN->giveNumber() == 0)
      lhs = this -> giveMassMatrix() -> GiveCopy() ;

   else {
      beta = domain -> giveTimeIntegrationScheme() -> giveBeta() ;
      if (beta == 0.0)
	 lhs = this -> giveMassMatrix() -> GiveCopy() ;
      else {
	 m   = this -> giveMassMatrix() ;
	 k   = this -> giveStiffnessMatrix() ;
	 dt  = stepN -> giveTimeIncrement() ;
	 lhs = k -> Times(beta*dt*dt) -> plus(m) ;}}

   return lhs ;
}


FloatArray*  Element :: ComputeNewmarkRhsAt (TimeStep* stepN)
   // Computes the contribution of the receiver to the right-hand side of the
   // linear system, using Newmark's formula.
{
   FloatMatrix *K ;
   FloatArray  *f,*d,*dPred,*rhs ;

   f = this -> ComputeLoadVectorAt(stepN) ;
   K = this -> giveStiffnessMatrix () ;

   if (stepN->giveNumber() == 0) {
      d     = this -> ComputeVectorOf ('d',stepN) ;
      rhs   = K -> Times(d->negated()) -> add(f) ;
      delete d ;}
   else {
      dPred = this -> ComputeVectorOf ('D',stepN) ;
      rhs   = K -> Times(dPred->negated()) -> add(f) ;
      delete dPred ;}

   delete f ;
   return rhs ;
}


int  Element :: computeNumberOfDofs ()
   // Returns the total number of dofs of the receiver's nodes.
{
   int i,n ;

   n = 0 ;
   for (i=1 ; i<=numberOfNodes ; i++)
      n += this -> giveNode(i) -> giveNumberOfDofs() ;
   return n ;
}


FloatArray*  Element :: ComputeResultingBodyForceAt (TimeStep* stepN)
   // Computes at stepN the resulting force due to all body loads that act
   // on the receiver. This force is used by the element for computing its
   // body load vector.
{
   int         i,n,nLoads ;
   BodyLoad*   load ;
   FloatArray  *force,*resultant ;

   resultant = new FloatArray(0) ;
   nLoads    = this -> giveBodyLoadArray() -> giveSize() ;
   for (i=1 ; i<=nLoads ; i++) {
      n     = bodyLoadArray -> at(i) ;
      load  = (BodyLoad*) domain->giveLoad(n) ;
      force = load -> ComputeForceOn(this,stepN) ;
      resultant -> add(force) ;
      delete force ;}

   if (resultant->giveSize() == 0) {
      delete resultant ;
      return NULL ;}
   else
      return resultant ;
}


FloatArray*  Element :: ComputeRhsAt (TimeStep* stepN)
   // Computes the contribution of the receiver to the right-hand side of the
   // linear system.
{
   TimeIntegrationScheme* scheme ;

   scheme = domain -> giveTimeIntegrationScheme() ;
   if (scheme -> isStatic())
      return  this -> ComputeStaticRhsAt (stepN) ;
   else if (scheme -> isNewmark())
      return  this -> ComputeNewmarkRhsAt(stepN) ;
   else {
      printf ("Error : unknown time integration scheme : %c\n",scheme) ;
      exit(0) ;}
   return new FloatArray(0);
}


FloatMatrix*  Element :: ComputeStaticLhsAt (TimeStep* stepN)
   // Computes the contribution of the receiver to the left-hand side of the
   // linear system, in a static analysis.
{
   return  this->giveStiffnessMatrix()->GiveCopy() ;
}


FloatArray*  Element :: ComputeStaticRhsAt (TimeStep* stepN)
   // Computes the contribution of the receiver to the right-hand side of the
   // linear system, in a static analysis.
{
   return  this->ComputeLoadVectorAt(stepN) ;
}


FloatMatrix*  Element :: computeStiffnessMatrix ()
   // Computes numerically the stiffness matrix of the receiver.
{
   int         i ;
   double      dV ;
   FloatMatrix *b,*db,*d ;
   GaussPoint  *gp ;

   stiffnessMatrix = new FloatMatrix() ;
   for (i=0 ; i<numberOfGaussPoints ; i++) {
      gp = gaussPointArray[i] ;
      b  = this -> ComputeBmatrixAt(gp) ;
      d  = this -> giveConstitutiveMatrix() ;
      dV = this -> computeVolumeAround(gp) ;
      db = d -> Times(b) ;
      stiffnessMatrix -> plusProduct(b,db,dV) ;

      delete b ;
      delete db ;}

   return stiffnessMatrix -> symmetrized() ;
}


FloatArray*  Element :: computeStrainVector (GaussPoint* gp, TimeStep* stepN)
   // Computes the vector containing the strains at the Gauss point gp of
   // the receiver, at time step stepN. The nature of these strains depends
   // on the element's type.
{
   FloatMatrix *b ;
   FloatArray  *u,*Epsilon ;

   b       = this -> ComputeBmatrixAt(gp) ;
   u       = this -> ComputeVectorOf('d',stepN) ;
   Epsilon = b -> Times(u) ;

   gp -> letStrainVectorBe(Epsilon) ;    // gp stores Epsilon, not a copy
   delete b ;
   delete u ;
   return Epsilon ;
}


FloatArray*  Element :: computeStressVector (GaussPoint* gp, TimeStep* stepN)
   // Computes the vector containing the stresses at the Gauss point gp of
   // the receiver, at time step stepN. The nature of these stresses depends
   // on the element's type.
{
   FloatArray *Epsilon,*Sigma ;

   Epsilon = gp -> giveStrainVector() ;
   if (! Epsilon)
      Epsilon = this->computeStrainVector(gp,stepN) ;
   Sigma = this->giveConstitutiveMatrix() -> Times(Epsilon) ;

   gp -> letStressVectorBe(Sigma) ;      // gp stores Sigma, not a copy
   return Sigma ;
}


FloatArray* Element :: ComputeVectorOf (char u, TimeStep* stepN)
   // Forms the vector containing the values of the unknown 'u' (e.g., the
   // displacement) of the dofs of the receiver's nodes.
{
   Node       *nodeI ;
   FloatArray *answer ;
   int        i,j,k,nDofs ;

   answer = new FloatArray(this->computeNumberOfDofs()) ;
   k      = 0 ;
   for (i=1 ; i<=numberOfNodes ; i++) {
      nodeI = this->giveNode(i) ;
      nDofs = nodeI->giveNumberOfDofs() ;
      for (j=1 ; j<=nDofs ; j++)
	 answer->at(++k) = nodeI->giveDof(j)->giveUnknown(u,stepN) ;}

   return answer ;
}


FloatArray* Element :: ComputeVectorOfPrescribed (char u, TimeStep* stepN)
   // Forms the vector containing the prescribed values of the unknown 'u'
   // (e.g., the prescribed displacement) of the dofs of the receiver's
   // nodes. Puts 0 at each free dof.
{
   Node       *nodeI ;
   Dof        *dofJ ;
   FloatArray *answer ;
   int        i,j,k,nDofs ;

   answer = new FloatArray(this->computeNumberOfDofs()) ;
   k      = 0 ;
   for (i=1 ; i<=numberOfNodes ; i++) {
      nodeI = this->giveNode(i) ;
      nDofs = nodeI->giveNumberOfDofs() ;
      for (j=1 ; j<=nDofs ; j++) {
	 dofJ = nodeI->giveDof(j) ;
	 if (dofJ -> hasBc())
	    answer->at(++k) = dofJ->giveUnknown(u,stepN) ;
	 else
	    answer->at(++k) = 0. ;}}

   return answer ;
}


IntArray*  Element :: giveBodyLoadArray ()
   // Returns the array which contains the number of every body load that act
   // on the receiver.
{
   int i,numberOfLoads ;

   if (! bodyLoadArray) {
      numberOfLoads = this -> readIfHas("bodyLoads") ;
      bodyLoadArray = new IntArray(numberOfLoads) ;
      for (i=1 ; i<=numberOfLoads ; i++)
	 bodyLoadArray->at(i) = this->readInteger("bodyLoads",i+1) ;}

   return bodyLoadArray ;
}


FloatMatrix*  Element :: giveConstitutiveMatrix ()
   // Returns the elasticity matrix {E} of the receiver.
{
   if (! constitutiveMatrix)
      this -> computeConstitutiveMatrix() ;

   return constitutiveMatrix ;
}


IntArray*  Element :: giveLocationArray ()
   // Returns the location array of the receiver. This array is obtained by
   // simply appending the location array of every node of the receiver.
{
   IntArray* nodalArray ;
   int       i ;

   if (! locationArray) {
      locationArray = new IntArray(0) ;
      for (i=1 ; i<=numberOfNodes ; i++) {
	 nodalArray    = this -> giveNode(i) -> giveLocationArray() ;
	 locationArray = locationArray -> followedBy(nodalArray) ;}}

   return locationArray ;
}


FloatMatrix*  Element :: giveMassMatrix ()
   // Returns the mass matrix of the receiver.
{
   if (! massMatrix)
      this -> computeMassMatrix() ;
   return massMatrix ;
}


Material*  Element :: giveMaterial ()
   // Returns the material of the receiver.
{
   if (! material)
      material = this -> readInteger("mat") ;

   return  domain -> giveMaterial(material) ;
}



Node*  Element :: giveNode (int i)
   // Returns the i-th node of the receiver.
{
   int  n ;

   if (! nodeArray)
      nodeArray = new IntArray(numberOfNodes) ;

   n = nodeArray->at(i) ;
   if (! n) {
      n = this -> readInteger("nodes",i) ;
      nodeArray->at(i) = n ;}

   return  domain -> giveNode(n) ;
}


FloatMatrix*  Element :: giveStiffnessMatrix ()
   // Returns the stiffness matrix of the receiver.
{
   if (! stiffnessMatrix)
      this -> computeStiffnessMatrix() ;
   return stiffnessMatrix ;
}


void  Element :: instanciateYourself ()
   // Gets from input file all data of the receiver.
{
   int i ;

#  ifdef VERBOSE
      printf ("instanciating element %d\n",number) ;
#  endif

   material = this -> readInteger("mat") ;
   nodeArray = new IntArray(numberOfNodes) ;
   for (i=1 ; i<=numberOfNodes ; i++)
      nodeArray->at(i) = this->readInteger("nodes",i) ;
   this -> giveBodyLoadArray() ;
}


Element*  Element :: ofType (char* aClass)
   // Returns a new element, which has the same number than the receiver,
   // but belongs to aClass (PlaneStrain, or Truss2D,..).
{
   Element* newElement ;

   if (! strncmp(aClass,"PN",2))
      newElement = new PlaneStrain(number,domain) ;
   else if (! strncmp(aClass,"T2",2))
      newElement = new Truss2D(number,domain) ;
   else {
      printf ("%s : unknown element type \n",aClass) ;
      exit(0) ;}

   return newElement ;
}


void  Element :: printOutputAt (TimeStep* stepN)
   // Performs end-of-step operations.
{
   int         i ;
   GaussPoint* gp ;
   FILE*       file ;

#  ifdef VERBOSE
      printf ("element %d printing output\n",number) ;
#  endif

   file = domain -> giveOutputStream() ;
   fprintf (file,"element %d :\n",number) ;

   for (i=1 ; i<=numberOfGaussPoints ; i++) {
      gp = gaussPointArray[i-1] ;
      this -> computeStrainVector(gp,stepN) ;
      this -> computeStressVector(gp,stepN) ;
      gp   -> printOutput() ;}
}


Element*  Element :: typed ()
   // Returns a new element, which has the same number than the receiver,
   // but is typed (PlaneStrain, or Truss2D,..).
{
   Element* newElement ;
   char     type[32] ;

   this -> readString("class",type) ;
   newElement = this -> ofType(type) ;

   return newElement ;
}


void  Element :: updateYourself ()
   // Updates the receiver at end of step.
{
   int i ;

#  ifdef VERBOSE
      printf ("updating element %d\n",number) ;
#  endif

   for (i=1 ; i<=numberOfGaussPoints ; i++)
      gaussPointArray[i-1] -> updateYourself() ;
}

