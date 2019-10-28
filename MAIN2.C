/*
	Program developed by Yves DUBOIS-PELERIN and Thomas ZIMMERMANN 
	at the Laboratory for Structural and Continuum mechanics, Ecole
	Polytechnique Federale de Lausanne under Swiss National	Science 
	Foundation grant number 21-27589.89, June 1992.
	
	Expurgated version developed by Yves Dubois-Pelerin for ZACE 
	Services Ltd, November 1992.
*/

//  MAIN
//  Solves finite element problems.
//
//  In order to measure accurately the time elapsed in the assembly phase
//  (formation of the linear system), this version of function 'main' clearly 
//  separates the input phase (instanciation) from the analysis phase.
//  Here the objects read anticipatively their data in the input file 
//  (using methods 'instanciateYourself()') ; by contrast, in the method 
//  'main' defined in file 'main.c', objects read their data only when 
//  required, so input and analysis phases overlap.
//  You may, or may not, use the compiling option DEBUG.

#include "domain.hxx"
#include "clock.h"


main ()
{
   Domain* mesh ;
   char    fileName[32] ;

   printf ("please enter the name of the data file : \n") ;
   gets (fileName) ;

   printf ("instanciation :            %.1f\n",timeNow()) ;
   mesh = new Domain(fileName) ;
   mesh -> instanciateYourself() ;
   printf ("start analysis :           %.1f\n",timeNow()) ;
   mesh -> solveYourself() ;
   delete mesh ;
   printf ("end of destruction :       %.1f\n",timeNow()) ;
}