//   **********************
//   *** CLASS MATERIAL ***
//   **********************
 

#include "femcmpnn.hxx"
#include "dictionr.hxx"


class Material : public FEMComponent
/*
   This class implements a material in a finite element problem. A material
   is an attribute of a domain. It is usually also attribute of many elements.
 DESCRIPTION
   The attribute 'propertyDictionary' contains all the properties of a mate-
   rial, like its Young modulus, its mass density, its area or its thickness.
 TASK
   Returning a material property (method 'give') ;
*/
{
   private :
      Dictionary*    propertyDictionary ;

   public :
      Material (int n,Domain* d) : FEMComponent(n,d)
				  { propertyDictionary = new Dictionary() ;}
      ~Material ()                { delete propertyDictionary ;}

      double   give (char) ;
      char*    giveClassName (char* s)       { return strcpy(s,"Material") ;}
      void     instanciateYourself () ;
      void     printYourself () ;
} ;









