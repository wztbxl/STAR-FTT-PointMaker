#ifndef STDBXMLWRITER_HH
#define STDBXMLWRITER_HH

#include "tableAcceptor.hh"
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <strstream.h>
#include <string.h>

//#include "TObject.h"

class StDbXmlWriter : public tableAcceptor {

protected:

  ostream* os;//!

public:

  StDbXmlWriter() {};
  StDbXmlWriter(ostream& ofs){ os=&ofs;};
  virtual ~StDbXmlWriter(){};

  virtual void streamHeader(const char* name);
  virtual void streamTableName(const char* name);
  virtual void streamEndTableName();
  virtual void streamAccessor();
  virtual void endAccessor();
  virtual void streamTail();

  virtual void pass(char* name, short& i, int len) ;  
  virtual void pass(char* name, int& i, int len);  
  virtual void pass(char* name, long& i, int len);  
  virtual void pass(char* name, unsigned short& i, int len) ;  
  virtual void pass(char* name, unsigned int& i, int len) ;  
  virtual void pass(char* name, unsigned long& i, int len) ;  

  virtual void pass(char* name, float& i, int len);
  virtual void pass(char* name, double& i, int len);
  virtual void pass(char* name, char* i, int len);
  virtual void pass(char* name, unsigned char* i, int len) ;
  virtual void pass(char* name, short* i, int len) ;  
  virtual void pass(char* name, int* i, int len);  
  virtual void pass(char* name, long* i, int len);  
  virtual void pass(char* name, unsigned short* i, int len) ;  
  virtual void pass(char* name, unsigned int* i, int len) ;  
  virtual void pass(char* name, unsigned long* i, int len) ;  
  virtual void pass(char* name, float* i, int len);
  virtual void pass(char* name, double* i, int len);

  //ClassDef(StDbXmlWriter,1)

};

#endif


