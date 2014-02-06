#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "exception.h"

Exception::Exception()
{
  type = UnknownException;
  sprintf(message,"An unknown error occured.");
}

Exception::Exception(const char *message, ...)
{
  type = UnknownException;
  this->message = new char[512];

  va_list arg;
  va_start(arg, message);
  vsprintf(this->message, message, arg);
  va_end(arg);
}

Exception::~Exception()
{
  delete [] message;
}

void Exception::Report()
{
	cout << "Error! " << message << endl;
}


FileNotFoundException::FileNotFoundException()
{
  type = FileNotFound;
  sprintf(message,"An unknown file could not be found.");
}

FileNotFoundException::FileNotFoundException(char *fname)
{
  type = FileNotFound;
  this->message = new char[512];
  sprintf(message, "The file '%s' could not be found.", fname);
}


ObjectNotFoundException::ObjectNotFoundException()
{
  type = ObjectNotFound;
  sprintf(message, "An unknown object could not be retrieved.");
}

ObjectNotFoundException::ObjectNotFoundException(const char *message, ...)
{
  type = ObjectNotFound;
  this->message = new char[512];

  va_list arg;
  va_start(arg, message);
  vsprintf(this->message, message, arg);
  va_end(arg);
}
