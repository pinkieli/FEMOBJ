//   file TRUSS2D.CXX

#include "truss2d.hxx"
#include "domain.hxx"
#include "node.hxx"
#include "material.hxx"
#include "gausspnt.hxx"
#include "flotmtrx.hxx"
#include "diagmtrx.hxx"
#include "flotarry.hxx"
#include "intarray.hxx"
#include <stdlib.h>
#include <math.h>


Truss2D :: Truss2D (int n, Domain* aDomain)
	 : Element (n,aDomain)
   // Constructor.
{
   numberOfNodes       = 2 ;
   rotationMatrix      = NULL ;
   length              = 0. ;
   pitch               = 10. ;			// a dummy value
   this -> computeGaussPoints() ;
}


FloatMatrix*  Truss2D :: ComputeBmatrixAt (GaussPoint* aGaussPoint)
   // Returns the matrix  { -1/L  0  1/L  0 } .
{
   double       coeff ;
   FloatMatrix* answer ;

   coeff  = 1. / this->giveLength() ;
   answer = new FloatMatrix(1,4) ;
   answer->at(1,1) = 0. - coeff ;
   answer->at(1,2) = 0. ;
   answer->at(1,3) = coeff ;
   answer->at(1,4) = 0. ;

   return answer ;
}


FloatArray*  Truss2D :: ComputeBodyLoadVectorAt (TimeStep* stepN)
   // Returns the body load vector of the receiver, at stepN, expressed in
   // the global axes.
{
   FloatArray*  b ;
   FloatMatrix* R ;

   b = this -> Element::ComputeBodyLoadVectorAt(stepN) ;
   if (b) {
      R = this -> giveRotationMatrix() ;
      return  b -> rotatedWith(R,'t') ;}
   else
      return NULL ;
}


FloatMatrix*  Truss2D :: computeConstitutiveMatrix ()
   // Returns a pointer to the 1-by-1 elasticity matrix {E} of the receiver.
{
   Material* mat ;
   double    e,a ;

   mat = this->giveMaterial() ;
   e   = mat->give('E') ;
   a   = mat->give('A') ;
   constitutiveMatrix          = new FloatMatrix(1,1) ;
   constitutiveMatrix->at(1,1) = e*a ;
   return constitutiveMatrix ;
}


void  Truss2D :: computeGaussPoints ()
   // Sets up the array of Gauss Points of the receiver.
{
   FloatArray* coord ;

   numberOfGaussPoints = 1 ;
   gaussPointArray     = new GaussPoint* [numberOfGaussPoints] ;

   coord              = new FloatArray(1) ;
   coord->at(1)       = 0. ;
   gaussPointArray[0] = new GaussPoint(this,1,coord,2.) ;
}


FloatMatrix*  Truss2D :: computeMassMatrix ()
   // Returns the lumped mass matrix of the receiver. This expression is
   // valid in both local and global axes.
{
   Material* mat ;
   double    halfMass ;

   mat        = this -> giveMaterial() ;
   halfMass   = mat->give('d') * mat->give('A') * this->giveLength() / 2.;
   massMatrix = new DiagonalMatrix(4) ;
   massMatrix -> at(1,1) = halfMass ;
   massMatrix -> at(2,2) = halfMass ;
   massMatrix -> at(3,3) = halfMass ;
   massMatrix -> at(4,4) = halfMass ;

   return massMatrix ;
}


FloatMatrix*  Truss2D :: ComputeNmatrixAt (GaussPoint* aGaussPoint)
   // Returns the displacement interpolation matrix {N} of the receiver, eva-
   // luated at aGaussPoint.
{
   double       ksi,n1,n2 ;
   FloatMatrix* answer ;

   ksi = aGaussPoint -> giveCoordinate(1) ;
   n1  = (1. - ksi) * 0.5 ;
   n2  = (1. + ksi) * 0.5 ;
   answer = new FloatMatrix(2,4) ;

   answer->at(1,1) = n1 ;
   answer->at(1,3) = n2 ;
   answer->at(2,2) = n1 ;
   answer->at(2,4) = n2 ;

   return answer ;
}


FloatArray*  Truss2D :: ComputeResultingBodyForceAt (TimeStep* stepN)
{
   FloatMatrix  *rot ;
   FloatArray   *f ;
   double       cos,sin,fX,fY,area ;

   f   = this -> Element::ComputeResultingBodyForceAt(stepN) ;
   if (f) {
      rot = this -> giveRotationMatrix() ;
      cos = rot->at(1,1) ;
      sin = rot->at(1,2) ;
      fX  = f->at(1) ;
      fY  = f->at(2) ;
      f->at(1) =  cos*fX + sin*fY ;
      f->at(2) = -sin*fX + cos*fY ;
      area = this->giveMaterial()->give('A') ;
      return f->times(area) ;}
   else
      return NULL ;
}


FloatMatrix*  Truss2D :: computeStiffnessMatrix ()
   // Returns the stiffness matrix of the receiver, expressed in the global
   // axes.
{
   Material* mat ;
   double    coeff ;

   mat   = this -> giveMaterial() ;
   coeff = mat->give('E') * mat->give('A') / this->giveLength() ;
   stiffnessMatrix = new FloatMatrix(4,4) ;
   stiffnessMatrix->at(1,1) =  coeff ;
   stiffnessMatrix->at(1,3) = -coeff ;
   stiffnessMatrix->at(3,1) = -coeff ;
   stiffnessMatrix->at(3,3) =  coeff ;

   this -> giveRotationMatrix() ;
   return  stiffnessMatrix -> rotatedWith(rotationMatrix) ;
}


FloatArray*  Truss2D :: computeStrainVector (GaussPoint* gp, TimeStep* stepN)
   // Returns the vector containing the strain of the receiver. The method
   // inherited from Element cannot fit : an axis rotation is required here.
{
   FloatMatrix *b,*rot ;
   FloatArray  *u,*epsilon ;

   b       = this -> ComputeBmatrixAt (gp) ;
   u       = this -> ComputeVectorOf ('d',stepN) ;
   rot     = this -> giveRotationMatrix () ;
   epsilon = b -> Times (u->rotatedWith(rot,'n')) ;

   gp -> letStrainVectorBe(epsilon) ;       // gp stores epsilon, not a copy
   delete b ;
   delete u ;
   return epsilon ;
}


double  Truss2D :: computeVolumeAround (GaussPoint* aGaussPoint)
   // Returns the length of the receiver. This method is valid only if 1
   // Gauss point is used.
{
   return this->giveLength() ;
}


double  Truss2D :: giveLength ()
   // Returns the length of the receiver.
{
   double dx,dy ;
   Node   *nodeA,*nodeB ;

   if (length == 0.) {
      nodeA   = this->giveNode(1) ;
      nodeB   = this->giveNode(2) ;
      dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1) ;
      dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2) ;
      length  = sqrt(dx*dx + dy*dy) ;}

   return length ;
}


double  Truss2D :: givePitch ()
   // Returns the pitch of the receiver.
{
   double xA,xB,yA,yB ;
   Node   *nodeA,*nodeB ;

   if (pitch == 10.) {                // 10. : dummy initialization value
      nodeA  = this -> giveNode(1) ;
      nodeB  = this -> giveNode(2) ;
      xA     = nodeA->giveCoordinate(1) ;
      xB     = nodeB->giveCoordinate(1) ;
      yA     = nodeA->giveCoordinate(2) ;
      yB     = nodeB->giveCoordinate(2) ;
      pitch  = atan2(yB-yA,xB-xA) ;}

   return pitch ;
}


FloatMatrix*  Truss2D :: giveRotationMatrix ()
   // Returns the rotation matrix of the receiver.
{
   double sine,cosine ;

   if (! rotationMatrix) {
      sine           = sin (this->givePitch()) ;
      cosine         = cos (pitch) ;
      rotationMatrix = new FloatMatrix(4,4) ;
      rotationMatrix -> at(1,1) =  cosine ;
      rotationMatrix -> at(1,2) =  sine   ;
      rotationMatrix -> at(2,1) = -sine   ;
      rotationMatrix -> at(2,2) =  cosine ;
      rotationMatrix -> at(3,3) =  cosine ;
      rotationMatrix -> at(3,4) =  sine   ;
      rotationMatrix -> at(4,3) = -sine   ;
      rotationMatrix -> at(4,4) =  cosine ;}

   return rotationMatrix ;
}
