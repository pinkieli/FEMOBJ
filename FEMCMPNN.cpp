//   file FEMCMPNN.CXX
 
#include "femcmpnn.hxx"
#include "domain.hxx"
#include "string.hxx"
#include "freader.hxx"                /* SIZE is defined in this file */
#include <stdio.h>
#include <stdlib.h>


double  FEMComponent :: read (char c)
   // Reads in the data file the value of the data 'c'.
{
   char buffer[2] ;

   buffer[0] = c ;
   buffer[1] = '\0' ;

   return  this -> read(buffer) ;
}


double  FEMComponent :: read (char* data, int i)
   // Reads in the data file the i-th value of 'data'.
{
   char n[8],value[32],className[32] ;

   this -> giveClassName(className) ;
   itos (number,n) ;
   domain -> giveInputStream() -> read(className,n,data,i,value) ;

   if (isNumber(value))
      return atof(value) ;
   else {
      printf ("Error when reading the value No %d of data '%s' of %s %d\n",
	       i,data,className,number) ;
      exit(0) ;}
   return 0.;
}


int  FEMComponent :: readIfHas (char* data)
   // If the receiver has 'data' in the data file, returns the (integral)
   // value of 'data', else returns 0.
{
   char n[8],value[32],className[32] ;
   int  answer ;                

   this->giveClassName(className) ;
   itos (number,n) ;
   answer = domain->giveInputStream()->readIfHas(className,n,data,value);
   if (answer)
      answer = atoi(value) ;

   return answer ;
}


int  FEMComponent :: readInteger (char* data, int i)
   // Reads in the data file the i-th value of 'data'.
{
   char n[8],value[32],className[32] ;

   this -> giveClassName(className) ;
   itos (number,n) ;
   domain -> giveInputStream() -> read(className,n,data,i,value) ;

   if (isInteger(value))
      return atoi(value) ;
   else {
      printf ("Error when reading the value No %d of data '%s' of %s %d\n",
	       i,data,className,number) ;
      exit(0) ;}
   return 0;
}


int  FEMComponent :: readNumberOf (char* k)
   // Reads in the data file the value following the keyword 'k'.
{
   char value[32] ;

   domain -> giveInputStream() -> read(k,value) ;
   return atoi(value) ;
}


void  FEMComponent :: readString (char* data, int i, char* s)
   // Reads in the data file the i-th value of 'data'. Initializes 's' to it.
{
   char n[8],className[32] ;

   this -> giveClassName(className) ;
   itos (number,n) ;
   domain -> giveInputStream() -> read(className,n,data,i,s) ;
}


int  FEMComponent :: readWhetherHas (char* data)
   // If the receiver has the real (ie floating point) data 'data' in the
   // data file, returns True, else returns False.
   // This method is used by objects when attempting to read their attributes
   // from the data file. It does not have the same aspect as readIfHas(),
   // since the value of a real data may be 0.0, undistinguishable from False.
{
   char n[8],value[32],className[32] ;
   int  answer ;

   this -> giveClassName(className) ;
   itos (number,n) ;
   answer = domain->giveInputStream()->readIfHas(className,n,data,value) ;

   return answer ;
}

