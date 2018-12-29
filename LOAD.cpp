//   file LOAD.CXX

#include "load.hxx"
#include "deadwght.hxx"
#include "nodload.hxx"
#include "boundary.hxx"
#include "initial.hxx"
#include "debug.def"
#include "verbose.def"
#include <stdlib.h>


Load :: Load (int i, Domain* aDomain)
	   : FEMComponent (i,aDomain)
   // Constructor. Creates a load with number i, belonging to aDomain.
{
   componentArray   = NULL ;
   loadTimeFunction = 0 ;
}


FloatArray*  Load :: giveComponentArray ()
   // Returns the array that contains the components of the receivert. If
   // this array does not exist yet, forms it by reading its values in the
   // data file.
{
   int i,size ;

   if (! componentArray) {
      size           = this -> readInteger("components") ;
      componentArray = new FloatArray(size) ;
      for  (i=1 ; i<=size ; i++)
	 componentArray->at(i) = this->read("components",i+1) ;}

   return componentArray ;
}


LoadTimeFunction*  Load :: giveLoadTimeFunction ()
   // Returns the load-time function of the receiver. Reads its number in the
   // data file if has not been done yet.
{
   if (! loadTimeFunction)
      loadTimeFunction = this -> readInteger("loadTimeFunction") ;

   return  domain -> giveLoadTimeFunction(loadTimeFunction) ;
}


void  Load :: instanciateYourself ()
{
#  ifdef VERBOSE
      printf ("instanciating load %d\n",number) ;
#  endif

   loadTimeFunction = this->readInteger("loadTimeFunction") ;
   this -> giveComponentArray() ;
}


Load*  Load :: ofType (char* aClass)
   // Returns a new load, which has the same number than the receiver,
   // but belongs to aClass (NodalLoad, DeadWeight,..).
{
   Load* newLoad ;

   if (! strncmp(aClass,"BoundaryCondition",5))
      newLoad = new BoundaryCondition(number,domain) ;
   else if (! strncmp(aClass,"DeadWeight",5))
      newLoad = new DeadWeight(number,domain) ;
   else if (! strncmp(aClass,"InitialCondition",5))
      newLoad = new InitialCondition(number,domain) ;
   else if (! strncmp(aClass,"NodalLoad",5))
      newLoad = new NodalLoad(number,domain) ;
   else {
      printf ("%s : unknown load type \n",aClass) ;
      exit(0) ;}

   return newLoad ;
}


Load*  Load :: typed ()
   // Returns a new load, which has the same number than the receiver,
   // but is typed (BoundaryCondition,DeadWeight,..).
{
   Load*  newLoad ;
   char   type[32] ;

   this -> readString("class",type) ;
   newLoad = this -> ofType(type) ;

   return newLoad ;
}











