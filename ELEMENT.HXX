//   *********************
//   *** CLASS ELEMENT *** 
//   *********************
 

#ifndef element_hxx

#include "femcmpnn.hxx"
#include "domain.hxx"
#include "flotmtrx.hxx"
class TimeStep ; class Node ; class Material ; class GaussPoint ;
class FloatMatrix ; class FloatArray ; class IntArray ;


class Element : public FEMComponent
/*
   This abstract class is the most important class of the program. It is the
   superclass of all classes implementing finite elements (bar, shell, etc).
   An element is an attribute of a domain.
 DESCRIPTION :
   The basic data of an element are the numbers of its 'numberOfNodes' nodes,
   stored in 'nodeArray', of its 'material', of its body loads (eg, the dead
   weight) stored in 'bodyLoadArray'. These data are obtained from the domain.
   The element also possesses 'numberOfGaussPoints' Gauss points, stored in
   'gaussPointArray'.
   The calculated data of an element are its 'massMatrix', its 'stiffnessMa-
   trix', its 'locationArray'. Since the load vector is recalculated at every
   time step, it is not given the status of attribute.
 TASKS :
   -defining itself :
     .typing itself (methods 'typed' and 'ofType'). When the domain creates
      an element, it actually creates a temporary instance of Element, then
      asks this element to transform itself into an element of the right type
      (PlaneStrain, Truss2D, etc) ;
     .obtaining its basic data : the element reads in the data file the num-
      ber of these objects, then obtains the data from the domain (methods
      'giveNode', 'giveMaterial',etc) ;
   -calculating its contribution to the problem :
     .calculating its mass matrix M, its stiffness matrix K, its load vector
      f, its location array ;
     .calculating its contribution to the LHS and RHS of the linear system,
      using Static,Newmark,etc, formula. These contributions are usually
      combinations of M,K,f.
     .requesting from the LHS and the RHS of the linear system to assemble
      these contributions, by the means of the location array ;
   -performing end-of-step operations :
     .calculating the strains and stresses at its Gauss points ;
     .printing its output in the data file and updating itself ;
*/
{
   protected :
      int           numberOfNodes ;
      IntArray*     nodeArray ;
      int           material ;
      IntArray*     bodyLoadArray ;
      int           numberOfGaussPoints ;
      GaussPoint**  gaussPointArray ;
      IntArray*     locationArray ;
      FloatMatrix*  massMatrix ;
      FloatMatrix*  stiffnessMatrix ;
      FloatMatrix*  constitutiveMatrix ;

   public :
      Element (int,Domain*) ;              // constructor
      virtual ~Element () ;                // destructor

      // assembly
      void                  assembleYourselfAt (TimeStep*) ;
      void                  assembleLhsAt (TimeStep*) ;
      void                  assembleRhsAt (TimeStep*) ;
      IntArray*             giveLocationArray () ;

      // left-hand side
      FloatMatrix*          ComputeLhsAt (TimeStep*) ;
      FloatMatrix*          ComputeNewmarkLhsAt (TimeStep*) ;
      FloatMatrix*          ComputeStaticLhsAt (TimeStep*) ;

      // right-hand side
      FloatArray*           ComputeRhsAt (TimeStep*) ;
      FloatArray*           ComputeNewmarkRhsAt (TimeStep*) ;
      FloatArray*           ComputeStaticRhsAt (TimeStep* stepN) ;

      // mass and stiffness matrix
      FloatMatrix*          giveMassMatrix () ;
      FloatMatrix*          giveStiffnessMatrix () ;
      FloatMatrix*          giveConstitutiveMatrix () ;
      virtual FloatMatrix*  computeMassMatrix () ;
      virtual FloatMatrix*  ComputeConsistentMassMatrix () ;
      virtual FloatMatrix*  computeStiffnessMatrix () ;
      virtual FloatMatrix*  computeConstitutiveMatrix ()    { return NULL ;}

      // load vector
      FloatArray*           ComputeLoadVectorAt (TimeStep*) ;
      FloatArray*           ComputeBcLoadVectorAt (TimeStep*) ;
      virtual FloatArray*   ComputeBodyLoadVectorAt (TimeStep*) ;
      virtual FloatArray*   ComputeResultingBodyForceAt (TimeStep*) ;

      // strains and stresses
      virtual FloatArray*   computeStrainVector (GaussPoint*,TimeStep*) ;
      FloatArray*           computeStressVector (GaussPoint*,TimeStep*) ;

      // vector of nodal unknowns
      virtual FloatArray*   ComputeVectorOf (char,TimeStep*) ;
      FloatArray*           ComputeVectorOfPrescribed (char,TimeStep*) ;
      int                   computeNumberOfDofs () ;

      // interpolation, numerical integration
      virtual void          computeGaussPoints ()             {}
      virtual FloatMatrix*  ComputeBmatrixAt (GaussPoint*)    {return NULL ;}
      virtual FloatMatrix*  ComputeNmatrixAt (GaussPoint*)    {return NULL ;}
      virtual double        computeVolumeAround (GaussPoint*) {return NULL ;}

      // data management
      Node*                 giveNode (int) ;
      Material*             giveMaterial () ;
      IntArray*             giveBodyLoadArray () ;
      void                  instanciateYourself () ;

      // time step termination
      void                  printOutputAt (TimeStep*) ;
      void                  updateYourself () ;

      // definition
      Element*              typed () ;
      Element*              ofType (char*) ;
      char*                 giveClassName (char* s)
                                      { return strcpy(s,"Element") ;}
      int                   giveNumber ()
				      { return FEMComponent::giveNumber() ;}
} ;

#define element_hxx
#endif








