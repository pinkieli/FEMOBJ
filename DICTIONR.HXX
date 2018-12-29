//   ************************
//   *** CLASS DICTIONARY ***
//   ************************
 

#ifndef dictionr_hxx

#include "pair.hxx"
#include <stdio.h>


class Dictionary
/*
   This class implements a linked list whose entries are Pairs (see below).
   Dictionaries are typically used by degrees of freedom for storing their
   unknowns.
 DESCRIPTION :
   A dictionary stores its pairs in a linked list form. It knows the first
   pair (attribute 'first') of the list. It also knows the last one (attri-
   bute 'last') in order to append fastly an additional pair.
 TASK :
   Storing pairs (method 'add') and returning the value of a pair (method
   'at').
*/
{
   enum { FALSE , TRUE } ;

   protected :
      Pair*  first ;
      Pair*  last ;

   public :
      Dictionary ()  { first=NULL ; last=NULL ;}
      ~Dictionary () ;

      Pair*        add (char,double) ;
      double&      at (char) ;
      int          includes (char) ;
      void         printYourself () ;
} ;

#define dictionr_hxx
#endif








