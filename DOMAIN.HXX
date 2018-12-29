//   ********************
//   *** CLASS DOMAIN ***
//   ********************

 
#ifndef domain_hxx

#include "list.hxx"
#include <stdio.h>
class Element ; class Node ; class Material ; class TimeIntegrationScheme ;
class TimeStep ; class Load ; class LoadTimeFunction ; class LinearSystem ;
class FileReader ;


class Domain
/*
   This class implements the domain of a finite element problem.
 DESCRIPTION :
   The domain stores its components (elements, nodes, materials, loads, load-
   time functions) in lists : 'elementList', 'nodeList', etc. Its component
   'timeIntegrationScheme' needs no list, since it is unique. The domain
   also possesses the system of linear algebraic equations of the problem
   ('linearSystem').
   The domain possesses two streams 'inputStream' and 'outputStream' that
   the components use in order to read their input and write their output in
   the data file.
 TASKS :
   - User interfacing. The domain is the primary interlocutor of the user.
     In particular, it understands the message 'solveYourself', which trig-
     gers the solving process.
   - Managing its components and the linear system.
     The domain maintains lists of its components. The domain is responsible
     for creating these components (methods 'giveElement','giveLinearSystem',
     etc). Asking the domain is the only way for a component (e.g., an ele-
     ment) to have access to another component (e.g., its material), or to
     the data file, or to the linear system.
   - Returning the input/output streams. Since these two streams act on the
     same file (the data file), the domain is responsible for always closing
     one when the other one is to be used.
*/
{
   private :
      char*       	      dataFileName ;
      List*        	      elementList ;
      List*        	      nodeList ;
      List*        	      materialList ;
      List*        	      loadList ;
      List*                   loadTimeFunctionList ;
      TimeIntegrationScheme*  timeIntegrationScheme ;
      int                     numberOfElements ;
      LinearSystem*  	      linearSystem ;
      FileReader*  	      inputStream ;
      FILE*      	      outputStream ;

   public :
      Domain () ;                             // constructors
      Domain (char*) ;
      ~Domain () ;                            // destructor

      // solving
      void               solveYourself () ;
      void               solveYourselfAt (TimeStep*) ;
      void               formTheSystemAt (TimeStep*) ;
      int                giveNumberOfElements () ;
      void               terminate (TimeStep*) ;

      // management of the mesh components
      Element*           giveElement (int) ;
      LinearSystem*      giveLinearSystem ()     { return linearSystem ;}
      Load*              giveLoad (int) ;
      LoadTimeFunction*  giveLoadTimeFunction (int) ;
      Material*          giveMaterial (int) ;
      Node*              giveNode (int) ;
      TimeIntegrationScheme*  giveTimeIntegrationScheme () ;
      void               instanciateYourself () ;

      // input / output
      char*              giveDataFileName () ;
      FileReader*        giveInputStream () ;
      FILE*              giveOutputStream () ;
      int                readNumberOf (char*) ;
} ;

#define domain_hxx
#endif









