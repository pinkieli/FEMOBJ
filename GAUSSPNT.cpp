//   file GAUSSPNT.CXX
 
#include "gausspnt.hxx"
#include "element.hxx"
#include "domain.hxx"
#include "string.hxx"
#include "debug.def"


GaussPoint :: GaussPoint (Element* e, int n, FloatArray* a, double w)
   // Constructor. Creates a Gauss point belonging to element e, with number
   // n, with coordinates a, with weight w.
{
   element      = e ;
   number       = n ;
   coordinates  = a ;
   weight       = w ;
   strainVector = NULL ;
   stressVector = NULL ;
}


GaussPoint :: ~GaussPoint ()
   // Destructor.
{
   delete coordinates ;
   delete strainVector ;
   delete stressVector ;
}


void  GaussPoint :: printOutput ()
   // Prints the strains and stresses on the data file.
{
   FILE *File ;
   int  i,n ;

   File = element -> giveDomain() -> giveOutputStream() ;

   fprintf (File,"  GP %d :  strains ",number) ;
   n = strainVector -> giveSize() ;
   for (i=1 ; i<=n ; i++)
      fprintf (File," % .4e",strainVector->at(i)) ;

   fprintf (File,"\n          stresses",number) ;
   n = stressVector -> giveSize() ;
   for (i=1 ; i<=n ; i++)
      fprintf (File," % .4e",stressVector->at(i)) ;

   fprintf (File,"\n") ;
}


void  GaussPoint :: updateYourself ()
   // Performs end-of-step updates.
{
   delete strainVector ;
   delete stressVector ;

   strainVector = NULL ;
   stressVector = NULL ;
}
 
