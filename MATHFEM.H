// This file contains some functions used in the finite element
// application.
// ref : Lippman p 104
// ! friend functions may not be declared inline (C++\Views p 17)
 

inline  int min (int i,int j)
   { return (i<=j ? i : j) ; } 

inline  int max (int i,int j)
   { return (i>=j ? i : j) ; }









