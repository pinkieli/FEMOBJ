//   file MATERIAL.CXX
 
#include "material.hxx"
#include "domain.hxx"
#include "verbose.def"
#include <stdlib.h>


double  Material :: give (char aProperty)
   // Returns the value of the property aProperty (e.g. the Young's modulus
   // 'E') of the receiver.
{
   double  value ;
 
   if (propertyDictionary -> includes(aProperty))
      value = propertyDictionary -> at(aProperty) ;
   else {                                         // read and store the value
      value = this -> read(aProperty) ;
      propertyDictionary -> add(aProperty,value) ;}
   return value ;
}


void  Material :: instanciateYourself ()
{
   double value ;

#  ifdef VERBOSE
      printf ("\ninstanciating material %d\n",number) ;
#  endif

   if (this->readWhetherHas("E")) {
      value = this -> read("E") ;
      propertyDictionary -> add('E',value) ;}

   if (this->readWhetherHas("n")) {
      value = this -> read("n") ;
      propertyDictionary -> add('n',value) ;}

   if (this->readWhetherHas("d")) {
      value = this -> read("d") ;
      propertyDictionary -> add('d',value) ;}

   if (this->readWhetherHas("A")) {
      value = this -> read("A") ;
      propertyDictionary -> add('A',value) ;}

   if (this->readWhetherHas("t")) {
      value = this -> read("t") ;
      propertyDictionary -> add('t',value) ;}
}


void  Material :: printYourself ()
   // Prints the receiver on screen.
{
   printf ("Material with properties : \n") ;
   propertyDictionary -> printYourself() ;
}


