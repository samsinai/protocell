#ifndef EXCEPTION
#define EXCEPTION

#include <iostream>
#include <string>

using namespace std;

/*
 * A very simple exception class that can be thrown uppon error.
 */
class Exception
{

 public:
  /*
   * Default constructor - unknown error.
   */
  Exception();

  /*
   * ...for use with formated message...
   */
  Exception(const char *message, ...);
 
  /*
   * Destructor; clean up message;
   */
  ~Exception();
  
  /*
   * Method called to display message in cout stream.
   */
  void Report();
  
  enum ExceptionType
  {
    UnknownException = 0,
    FileNotFound = 1,
    ObjectNotFound = 2
  };
  ExceptionType type;

 protected:
  /* the error message... */
  char* message;
};

class FileNotFoundException : public Exception
{
 public:
  FileNotFoundException();
  FileNotFoundException(char *fname);
};

class ObjectNotFoundException : public Exception 
{
 public:
  /*
   * Default constructor - unknown error.
   */
  ObjectNotFoundException();

  /*
   * ...for use with formated message...
   */
  ObjectNotFoundException(const char *message, ...);
};


#endif
