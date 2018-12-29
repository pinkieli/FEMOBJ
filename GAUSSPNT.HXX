//   *************************
//   *** CLASS GAUSS POINT ***
//   *************************
 

#include "flotarry.hxx"
class Element ;


class GaussPoint
/*
   This class implements a point for a gaussian quadrature. A Gauss point is
   usually attribute of an element.
 DESCRIPTION
   A Gauss point is identified by its 'number' and by 'element' - the element
   it belongs to. The array 'coordinates' defines its position in an axis sys-
   tem defined by its element. 'weight' is the pondaration of the point in the
   numerical integration ; 'weight' is naturally positive and smaller than 2.
   'strainVector' and 'stressVector' store the current values of strains and
   stresses ; the size of these arrays and the nature of their values is de-
   termined by the element.
 TASKS
   - returning its coordinates, its weight and the strains and stresses ;
   - reading/writing its attribute in a file.
 REMARK
   The Gauss point is a rather passive object : it does not compute its strains
   and stresses - it just stores them. They are computed by the element. Ifever
   materially nonlinearity is implemented, the role of the Gauss point is 
   likely to grow.
*/
{
   private :
      int          number ;
      Element*     element ;
      FloatArray*  coordinates ;
      double       weight ;
      FloatArray*  strainVector ;
      FloatArray*  stressVector ;

   public :
      GaussPoint  (Element*,int,FloatArray*,double) ;    // constructor
      ~GaussPoint  () ;                                  // destructor

      double       giveCoordinate (int i)      { return coordinates->at(i) ;}
      FloatArray*  giveCoordinates ()          { return coordinates ;}
      FloatArray*  giveStrainVector ()         { return strainVector ;}
      FloatArray*  giveStressVector ()         { return stressVector ;}
      double       giveWeight ()               { return weight ;}
      void         letStrainVectorBe (FloatArray* v) { strainVector = v ;}
      void         letStressVectorBe (FloatArray* v) { stressVector = v ;}
      void         printOutput () ;
      void         updateYourself () ;
} ;








