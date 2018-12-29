//   **********************	
//   *** CLASS TRUSS 2D ***
//   **********************


#include "element.hxx"


class Truss2D : public Element
/*
   This class implements a two-node truss bar element for two-dimensional
   analysis.
 DESCRIPTION :
   A truss bar element is characterized by its 'length' and its 'pitch'. The
   pitch is the angle in radians between the X-axis anf the axis of the
   element (oriented node1 to node2).
   The 'rotationMatrix' R is such that u{loc}=R*u{glob}.
 TASKS :
   - same as PlaneStrain element ;
   - expressing M,K,f,etc, in global axes. Methods like 'computeStiffness-
     Matrix' of class Element are here overloaded in order to account for
     rotational effects.
*/
{
   private :
      double        length ;
      double        pitch ;
      FloatMatrix*  rotationMatrix ;

   public :
      Truss2D (int,Domain*) ;                       // constructor
      ~Truss2D ()   { delete rotationMatrix ;}      // destructor

      FloatMatrix*  ComputeBmatrixAt (GaussPoint*) ;
      FloatArray*   ComputeBodyLoadVectorAt (TimeStep*) ;
      FloatMatrix*  computeConstitutiveMatrix () ;
      void          computeGaussPoints () ;
      FloatMatrix*  computeMassMatrix () ;
      FloatMatrix*  ComputeNmatrixAt (GaussPoint*) ;
      FloatArray*   ComputeResultingBodyForceAt (TimeStep*) ;
      FloatMatrix*  computeStiffnessMatrix () ;
      FloatArray*   computeStrainVector (GaussPoint*,TimeStep*) ;
      double        computeVolumeAround (GaussPoint*) ;
      double        giveLength () ;
      double        givePitch () ;
      FloatMatrix*  giveRotationMatrix () ;
} ;










