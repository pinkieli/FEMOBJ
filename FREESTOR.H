//   ********************************************
//   *** DYNAMIC MEMORY ALLOCATION PROCEDURES ***
//   ********************************************


#ifndef freestor_h


/*
   This file does not define a class. Rather, it provides a few procedures
   related with dynamic memory allocation.
*/

   double*  allocDouble (int) ;
   int*     allocInt (int) ;
   void     freeStoreError () ;
   void     freeInt (int*) ;
   void     freeDouble (double*) ;


#define freestor_h
#endif
