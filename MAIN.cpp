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
//  The DEBUG option MUST be used (check in file 'debug.def').
//  See also file 'main2.c'.

#include "domain.hxx"
#include "freestor.h"
#include "compiler.def"
//#include "stdafx.h"

int main ()
{
	/*
#ifdef BORLAND
   set_new_handler (&freeStoreError) ;     // prevents memory overflow
#endif
   */

   Domain* mesh ;

   mesh = new Domain() ;
   mesh -> solveYourself() ;
   delete mesh ;
   return 0;
}



