//   file DOMAIN.CXX
 
#include "domain.hxx"
#include "element.hxx"
#include "timestep.hxx"
#include "node.hxx"
#include "material.hxx"
#include "load.hxx"
#include "loadtime.hxx"
#include "timinteg.hxx"
#include "linsyst.hxx"
#include "string.hxx"
#include "freader.hxx"
#include "clock.h"
#include "debug.def"
#include "verbose.def"


Domain :: Domain ()
   // Constructor. Creates a new domain.
{
   dataFileName          = NULL ;
   elementList           = new List(0) ;
   nodeList              = new List(0) ;
   materialList          = new List(0) ;
   loadList              = new List(0) ;
   loadTimeFunctionList  = new List(0) ;
   timeIntegrationScheme = NULL ;
   inputStream  	 = NULL ;
   outputStream 	 = NULL ;
   numberOfElements      = 0 ;
   linearSystem          = new LinearSystem() ;
}


Domain :: Domain (char* s)
   // Constructor. Creates a new domain with data file 's'.
{
   dataFileName = new char[strlen(s)+1] ;
   strcpy (dataFileName,s) ;

   elementList           = new List(0) ;
   nodeList              = new List(0) ;
   materialList          = new List(0) ;
   loadList              = new List(0) ;
   loadTimeFunctionList  = new List(0) ;
   timeIntegrationScheme = NULL ;
   inputStream  	 = NULL ;
   outputStream 	 = NULL ;
   numberOfElements      = 0 ;
   linearSystem          = new LinearSystem() ;
}


Domain :: ~Domain ()
   // Destructor.
{
   delete dataFileName ;
   delete elementList ;
   delete nodeList ;
   delete materialList ;
   delete loadList ;
   delete loadTimeFunctionList ;
   delete timeIntegrationScheme ;
   delete linearSystem ;
   delete inputStream ;
   if (outputStream)
      fclose(outputStream) ;
}


void  Domain :: formTheSystemAt (TimeStep* stepN)
   // Assembles the system of linear equations, at the current time step.
{
   int i,nNodes ;

   this -> giveNumberOfElements() ;
   for (i=1 ; i<=numberOfElements ; i++)
      this -> giveElement(i) -> assembleYourselfAt(stepN) ;

   nNodes = this -> readNumberOf("Node") ;
   for (i=1 ; i<=nNodes ; i++)
      this -> giveNode(i) -> assembleYourLoadsAt(stepN) ;
}


char*  Domain :: giveDataFileName ()
   // Returns the name of the file containing the data of the problem.
{
   char s[64] ;

   if (! dataFileName) {
      printf ("please enter the name of the data file : \n") ;
      gets (s) ;
      dataFileName = new char[strlen(s)+1] ;
      strcpy (dataFileName,s) ;}

   return dataFileName ;
}


Element*  Domain :: giveElement (int n)
   // Returns the n-th element. Creates this element if it does not exist yet.
{
   Element* elem ;

   if (elementList -> includes(n))
      elem = (Element*) elementList->at(n) ;
   else {
      elem = (Element*) Element(n,this).typed() ;
      elementList -> put(n,elem) ;}

   return elem ;
}


FileReader*  Domain :: giveInputStream ()
   // Returns an input stream on the data file of the receiver.
{
   if (inputStream)
      return inputStream -> reset() ;

   else {
      if (outputStream) {              // flush output stream, if it exists
	 fclose (outputStream) ;
	 outputStream = NULL ;}
      inputStream = new FileReader(this->giveDataFileName()) ;
      return inputStream ;}
}


Load*  Domain :: giveLoad (int n)
   // Returns the n-th load. Creates this load if it does not exist yet.
{
   Load* load ;

   if (loadList -> includes(n))
      load = (Load*) loadList->at(n) ;
   else {
      load = (Load*) Load(n,this).typed() ;
      loadList -> put(n,load) ;}

   return load ;
}


LoadTimeFunction*  Domain :: giveLoadTimeFunction (int n)
   // Returns the n-th load-time function. Creates this fuction if it does
   // not exist yet.
{
   LoadTimeFunction* ltf ;

   if (loadTimeFunctionList -> includes(n))
      ltf = (LoadTimeFunction*) loadTimeFunctionList->at(n) ;
   else {
      ltf = (LoadTimeFunction*) LoadTimeFunction(n,this).typed() ;
      loadTimeFunctionList -> put(n,ltf) ;}

   return ltf ;
}


Material*  Domain :: giveMaterial (int n)
   // Returns the n-th material. Creates this material if it does not exist
   // yet.
{
   Material* mat ;

   if (materialList -> includes(n))
      mat = (Material*) materialList -> at(n) ;
   else {
      mat = new Material(n,this) ;
      materialList  -> put(n,mat) ;}

   return mat ;
}


Node*  Domain :: giveNode (int n)
   // Returns the n-th node. Creates this node if it does not exist yet.
{
   Node *node ;

   if (nodeList -> includes(n))
      node = (Node*) nodeList -> at(n) ;
   else {
      node = new Node(n,this) ;
      nodeList  -> put(n,node) ;}

   return node ;
}


int  Domain :: giveNumberOfElements ()
   // Returns the number of elements the problem consists of.
{
   if (! numberOfElements)
      numberOfElements = this -> readNumberOf("Element") ;
   return numberOfElements ;
}


FILE*  Domain :: giveOutputStream ()
   // Returns an output stream on the data file of the receiver.
{
   if (! outputStream) {
      if (inputStream) {                // flush input stream, if it exists
	 delete inputStream ;
	 inputStream = NULL ;}
      outputStream = fopen (dataFileName,"a") ;}

   return outputStream ;
}


TimeIntegrationScheme*  Domain :: giveTimeIntegrationScheme ()
   // Returns the time integration algorithm. Creates it if it does not
   // exist yet.
{
   TimeIntegrationScheme* scheme ;

   if (timeIntegrationScheme)
      return  timeIntegrationScheme ;
   else {
      scheme                = TimeIntegrationScheme(1,this).typed() ;
      timeIntegrationScheme = scheme ;
      return scheme ;}
}


void  Domain :: instanciateYourself ()
   // Creates all objects mentioned in the data file.
   // Exception : time step 2 and subsequent ones are not instanciated.
{
   int i,n ;

#  ifdef VERBOSE
      printf ("Reading all data from input file \n") ;
#  endif

   this -> giveTimeIntegrationScheme () ;

   n = this->readNumberOf("Node") ;
   nodeList -> growTo(n) ;
   for (i=1 ; i<=n ; i++)
      this -> giveNode(i) -> instanciateYourself() ;

   n = this->giveNumberOfElements() ;
   elementList -> growTo(n) ;
   for (i=1 ; i<=n ; i++)
      this -> giveElement(i) -> instanciateYourself() ;

   n = this->readNumberOf("Material") ;
   materialList -> growTo(n) ;
   for (i=1 ; i<=n ; i++)
      this -> giveMaterial(i) -> instanciateYourself() ;

   n = this->readNumberOf("Load") ;
   loadList -> growTo(n) ;
   for (i=1 ; i<=n ; i++)
      this -> giveLoad(i) -> instanciateYourself() ;

   n = this->readNumberOf("LoadTimeFunction") ;
   loadTimeFunctionList -> growTo(n) ;
   for (i=1 ; i<=n ; i++)
      this -> giveLoadTimeFunction(i) -> instanciateYourself() ;

   linearSystem -> carveYourselfFor(this) ;
}


int  Domain :: readNumberOf (char* type)
   // Gets from the data file the number of objects of type 'type' (e.g.
   // Element, or Node) that the receiver possesses.
{
   char value[8] ;

   this -> giveInputStream() -> read(type,value) ;
   return atoi(value) ;
}


void  Domain :: solveYourself ()
   // Solves the problem described by the receiver.
{
   TimeStep* currentStep ;

   this -> giveTimeIntegrationScheme() ;
   while (currentStep = timeIntegrationScheme->giveNextStep())
      this -> solveYourselfAt(currentStep) ;
}


void  Domain :: solveYourselfAt (TimeStep* stepN)
   // Solves the problem at the current time step.
{
   printf ("start forming the system : %.1f\n",timeNow()) ;
   this -> formTheSystemAt(stepN) ;

   printf ("start solving the system : %.1f\n",timeNow()) ;
     linearSystem -> solveYourself() ;

   printf ("start step termination   : %.1f\n",timeNow()) ;
     this -> terminate(stepN) ;

   printf ("end of step :              %.1f\n",timeNow()) ;
}


void  Domain :: terminate (TimeStep* stepN)
   // Performs all operations (printings and updates) for terminating time
   // step stepN.
{
   TimeStep* nextStep ;
   Element*  elem ;
   int       i,nNodes ;
   FILE*     File ;

   File = this -> giveOutputStream() ;
   fprintf (File,"\nSolution %d\n",stepN->giveNumber()) ;

   nNodes = nodeList -> giveSize() ;
   for (i=1 ; i<=nNodes ; i++)
      this -> giveNode(i) -> printOutputAt(stepN) ;

   for (i=1 ; i<=numberOfElements ; i++) {
      elem = this -> giveElement(i) ;
      elem -> printOutputAt(stepN) ;
      elem -> updateYourself() ;}

   nNodes = nodeList -> giveSize() ;
   for (i=1 ; i<=nNodes ; i++)
      this -> giveNode(i) -> updateYourself() ;

   if (stepN->isNotTheLastStep()) {
      nextStep = new TimeStep(stepN->giveNumber()+1,timeIntegrationScheme) ;
      if (nextStep->requiresNewLhs())
	 linearSystem -> updateYourself() ;
      else
	 linearSystem -> updateYourselfExceptLhs() ;
      delete nextStep ;}

   timeIntegrationScheme -> updateYourself() ;
}

