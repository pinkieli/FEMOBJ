// This file contains some functions used in the finite element
// application.
// ref : Lippman p 104
// ! friend functions may not be declared inline (C++\Views p 17)
 

inline  min (int i,int j)
   { return (i<=j ? i : j) ; } 

inline  max (int i,int j)
   { return (i>=j ? i : j) ; }









