//   **************************
//   *** CLASS FEMCOMPONENT ***
//   **************************

 
#ifndef femcmpnn_hxx

#include "string.hxx"
class Domain ; class FileReader ;


class FEMComponent
/*
   This class is an abstract class, the superclass of all classes that imple-
   ment the components of a finite element mesh : elements, nodes, time steps,
   materials, loads and load-time functions.
 DESCRIPTION
   This class defines the two attributes common to all component classes ;
   'number' is primarily used for reading data in the data file. 'domain' is
   used for communicating with other components (e.g., for an element to ob-
   tain its material), for accessing the linear system and the data file.
 TASKS
   Reading data in the data file.
*/
{
   protected :
      int      number ;
      Domain*  domain ;

   public :
      FEMComponent () {}                                    // constructors
      FEMComponent (int n,Domain* d) {number=n; domain=d;}
      virtual ~FEMComponent () {}                           // destructor

      virtual char*  giveClassName (char*) = 0 ;
      Domain*        giveDomain ()                 {return domain ;}
      int            giveNumber ()                 {return number ;}
      double         read (char) ;
      double         read (char* d)                {return this->read(d,1) ;}
      double         read (char*,int) ;
      int            readIfHas (char*) ;
      int            readInteger (char* d)    {return this->readInteger(d,1);}
      int            readInteger (char* d,int i);
      int            readNumberOf (char*) ;
      void           readString (char* d,char* s)  {this->readString(d,1,s) ;}
      void           readString (char*,int,char*);
      int            readWhetherHas (char*) ;
} ;

#define femcmpnn_hxx
#endif








