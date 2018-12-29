//   **************************
//   *** CLASS PLANE STRAIN ***
//   **************************
 

#include "element.hxx"
#include "polymtrx.hxx"


class PlaneStrain : public Element
/*
   This class implements an isoparametric four-node quadrilateral plane-
   strain elasticity finite element. Each node has 2 degrees of freedom.
 DESCRIPTION :
   One single additional attribute is needed for Gauss integration purpose :
   'jacobianMatrix'. This 2x2 matrix contains polynomials.
 TASKS :
   - calculating its Gauss points ;
   - calculating its B,D,N matrices and dV.
*/
{
   protected :
      PolynomialMatrix*  jacobianMatrix ;

   public :
      PlaneStrain (int,Domain*) ;                          // constructor
      ~PlaneStrain ()  {delete jacobianMatrix ;}           // destructor

      FloatMatrix*       ComputeBmatrixAt (GaussPoint*) ;
      FloatMatrix*       computeConstitutiveMatrix () ;
      void               computeGaussPoints () ;
      FloatMatrix*       ComputeNmatrixAt (GaussPoint*) ;
      double             computeVolumeAround (GaussPoint*) ;
      PolynomialMatrix*  giveJacobianMatrix () ;
} ;

   








