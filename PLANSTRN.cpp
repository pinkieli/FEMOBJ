//   file PLANSTRN.CXX
 
#include "planstrn.hxx"
#include "node.hxx"
#include "material.hxx"
#include "gausspnt.hxx"
#include "polynoxy.hxx"
#include "flotmtrx.hxx"
#include "diagmtrx.hxx"
#include "flotarry.hxx"
#include "intarray.hxx"
#include "domain.hxx"
#include <math.h>
#include <stdio.h>


PlaneStrain :: PlaneStrain (int n, Domain* aDomain)
	     : Element (n,aDomain)
   // Constructor.
{
   numberOfNodes  = 4 ;
   jacobianMatrix = NULL ;
   this -> computeGaussPoints() ;
}


FloatMatrix*  PlaneStrain :: ComputeBmatrixAt (GaussPoint *aGaussPoint)
   // Returns the [3x8] strain-displacement matrix {B} of the receiver, eva-
   // luated at aGaussPoint.
{
   FloatMatrix *jacGP,*inv,*answer ;
   FloatArray  *coord ;
   double      ksi,eta,ksiX,ksiY,etaX,etaY,
	       n1Ksi,n2Ksi,n3Ksi,n4Ksi,n1Eta,n2Eta,n3Eta,n4Eta ;

   // natural coordinates ksi and eta :
   coord = aGaussPoint -> giveCoordinates() ;
   ksi   = coord -> at(1) ;
   eta   = coord -> at(2) ;

   // partial derivatives of ksi and eta, with respect to x and y :
   jacGP = this -> giveJacobianMatrix() -> EvaluatedAt(coord) ;
   inv   = jacGP -> GiveInverse() ;

   ksiX = inv -> at(1,1) ;
   ksiY = inv -> at(1,2) ;
   etaX = inv -> at(2,1) ;
   etaY = inv -> at(2,2) ;

   // partial derivatives of the shape functions N_i, with respect to ksi
   // and eta :
   n1Ksi = (-1. + eta) * 0.25 ;
   n2Ksi = ( 1. - eta) * 0.25 ;
   n3Ksi = ( 1. + eta) * 0.25 ;
   n4Ksi = (-1. - eta) * 0.25 ;
   n1Eta = (-1. + ksi) * 0.25 ;
   n2Eta = (-1. - ksi) * 0.25 ;
   n3Eta = ( 1. + ksi) * 0.25 ;
   n4Eta = ( 1. - ksi) * 0.25 ;

   // B matrix  -  3 rows : epsilon-X, epsilon-Y, gamma-XY  :
   answer = new FloatMatrix(3,8) ;
   answer->at(1,1) = n1Ksi*ksiX + n1Eta*etaX ;
   answer->at(1,3) = n2Ksi*ksiX + n2Eta*etaX ;
   answer->at(1,5) = n3Ksi*ksiX + n3Eta*etaX ;
   answer->at(1,7) = n4Ksi*ksiX + n4Eta*etaX ;
   answer->at(2,2) = n1Ksi*ksiY + n1Eta*etaY ;
   answer->at(2,4) = n2Ksi*ksiY + n2Eta*etaY ;
   answer->at(2,6) = n3Ksi*ksiY + n3Eta*etaY ;
   answer->at(2,8) = n4Ksi*ksiY + n4Eta*etaY ;
   answer->at(3,1) = n1Ksi*ksiY + n1Eta*etaY ;
   answer->at(3,2) = n1Ksi*ksiX + n1Eta*etaX ;
   answer->at(3,3) = n2Ksi*ksiY + n2Eta*etaY ;
   answer->at(3,4) = n2Ksi*ksiX + n2Eta*etaX ;
   answer->at(3,5) = n3Ksi*ksiY + n3Eta*etaY ;
   answer->at(3,6) = n3Ksi*ksiX + n3Eta*etaX ;
   answer->at(3,7) = n4Ksi*ksiY + n4Eta*etaY ;
   answer->at(3,8) = n4Ksi*ksiX + n4Eta*etaX ;

   delete jacGP ;
   delete inv ;

   return answer ;
}

   
FloatMatrix*  PlaneStrain :: computeConstitutiveMatrix ()
   // Returns the [3x3] elasticity matrix {E} of the receiver.
{
   Material *mat ;
   double   e,nu,ee,shear ;

   mat   = this -> giveMaterial() ;
   e     = mat -> give('E') ;
   nu    = mat -> give('n') ;
   ee    = e / ((1.+nu) * (1.-nu-nu)) ;
   shear = e / (2.+nu+nu) ;
   constitutiveMatrix = new FloatMatrix(3,3) ;

   constitutiveMatrix->at(1,1) = (1.-nu) * ee ;
   constitutiveMatrix->at(1,2) =     nu  * ee ;
   constitutiveMatrix->at(2,1) =     nu  * ee ;
   constitutiveMatrix->at(2,2) = (1.-nu) * ee ;
   constitutiveMatrix->at(3,3) =  shear ;

   return constitutiveMatrix ;
}


void  PlaneStrain :: computeGaussPoints ()
   // Sets up the array containing the four Gauss points of the receiver.
{
   FloatArray *coord1,*coord2,*coord3,*coord4 ;
   double     a,aMinus ;

   a                   = 1. / sqrt(3.) ;
   aMinus              = -a ;
   numberOfGaussPoints = 4 ;
   gaussPointArray     = new GaussPoint* [numberOfGaussPoints] ;

   coord1             = new FloatArray(2) ;
   coord1 -> at(1)    = aMinus ;
   coord1 -> at(2)    = aMinus ;
   gaussPointArray[0] = new GaussPoint(this,1,coord1,1.) ;

   coord2             = new FloatArray(2) ;
   coord2 -> at(1)    = a ;
   coord2 -> at(2)    = aMinus ;
   gaussPointArray[1] = new GaussPoint(this,2,coord2,1.) ;

   coord3             = new FloatArray(2) ;
   coord3 -> at(1)    = a ;
   coord3 -> at(2)    = a ;
   gaussPointArray[2] = new GaussPoint(this,3,coord3,1.) ;

   coord4             = new FloatArray(2) ;
   coord4 -> at(1)    = aMinus ;
   coord4 -> at(2)    = a ;
   gaussPointArray[3] = new GaussPoint(this,4,coord4,1.) ;
}


FloatMatrix*  PlaneStrain :: ComputeNmatrixAt (GaussPoint* aGaussPoint)
   // Returns the displacement interpolation matrix {N} of the receiver, eva-
   // luated at aGaussPoint.
{
   double       ksi,eta,n1,n2,n3,n4 ;
   FloatMatrix* answer ;

   ksi = aGaussPoint -> giveCoordinate(1) ;
   eta = aGaussPoint -> giveCoordinate(2) ;
   n1  = (1. - ksi) * (1. - eta) * 0.25 ;
   n2  = (1. + ksi) * (1. - eta) * 0.25 ;
   n3  = (1. + ksi) * (1. + eta) * 0.25 ;
   n4  = (1. - ksi) * (1. + eta) * 0.25 ;
   answer = new FloatMatrix(2,8) ;

   answer->at(1,1) = n1 ;
   answer->at(1,3) = n2 ;
   answer->at(1,5) = n3 ;
   answer->at(1,7) = n4 ;
   answer->at(2,2) = n1 ;
   answer->at(2,4) = n2 ;
   answer->at(2,6) = n3 ;
   answer->at(2,8) = n4 ;

   return answer ;
}


double  PlaneStrain :: computeVolumeAround (GaussPoint* aGaussPoint)
   // Returns the portion of the receiver which is attached to aGaussPoint.
{
   FloatMatrix* jacob ;
   FloatArray*  coord ;
   double       determinant,weight,thickness,volume ;

   coord       = aGaussPoint -> giveCoordinates() ;
   jacob       = this -> giveJacobianMatrix() -> EvaluatedAt(coord) ;
   determinant = fabs (jacob->giveDeterminant()) ;
   weight      = aGaussPoint -> giveWeight() ;
   thickness   = this -> giveMaterial() -> give('t') ;

   volume      = determinant * weight * thickness ;

   delete jacob ;
   return volume ;
}


PolynomialMatrix*  PlaneStrain :: giveJacobianMatrix ()
   // Returns the jacobian matrix  J (x,y)/(ksi,eta)  of the receiver. Com-
   // putes it if it does not exist yet.
{
   Node         *node1,*node2,*node3,*node4 ;
   PolynomialXY *j11,*j12,*j21,*j22 ;
   double       x1,x2,x3,x4,y1,y2,y3,y4 ;

   if (! jacobianMatrix) {
      node1 = this -> giveNode(1) ;
      node2 = this -> giveNode(2) ;
      node3 = this -> giveNode(3) ;
      node4 = this -> giveNode(4) ;
      x1 = node1 -> giveCoordinate(1) ;
      x2 = node2 -> giveCoordinate(1) ;
      x3 = node3 -> giveCoordinate(1) ;
      x4 = node4 -> giveCoordinate(1) ;
      y1 = node1 -> giveCoordinate(2) ;
      y2 = node2 -> giveCoordinate(2) ;
      y3 = node3 -> giveCoordinate(2) ;
      y4 = node4 -> giveCoordinate(2) ;

      j11 = new PolynomialXY(1) ;
      j11 -> at(1) = (-x1 + x2 + x3 - x4) * 0.25 ;
      j11 -> at(2) =   0.                        ;
      j11 -> at(3) = ( x1 - x2 + x3 - x4) * 0.25 ;

      j12 = new PolynomialXY(1) ;
      j12 -> at(1) = (-x1 - x2 + x3 + x4) * 0.25 ;
      j12 -> at(2) = ( x1 - x2 + x3 - x4) * 0.25 ;
      j12 -> at(3) =   0.                        ;

      j21 = new PolynomialXY(1) ;
      j21 -> at(1) = (-y1 + y2 + y3 - y4) * 0.25 ;
      j21 -> at(2) =   0.                        ;
      j21 -> at(3) = ( y1 - y2 + y3 - y4) * 0.25 ;

      j22 = new PolynomialXY(1) ;
      j22 -> at(1) = (-y1 - y2 + y3 + y4) * 0.25 ;
      j22 -> at(2) = ( y1 - y2 + y3 - y4) * 0.25 ;
      j22 -> at(3) =   0.                        ;

      jacobianMatrix = new PolynomialMatrix(2,2) ;
      jacobianMatrix -> at(1,1) = j11 ;
      jacobianMatrix -> at(1,2) = j12 ;
      jacobianMatrix -> at(2,1) = j21 ;
      jacobianMatrix -> at(2,2) = j22 ;}

   return jacobianMatrix ;
}


