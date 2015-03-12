#include <iostream>
#include <fstream>
#include <fcntl.h>
#include "exception.h"
#include "fileioutils.h"

using namespace std;

void FileIOUtils::NotifyFileReading(char* fname)
{
  cout << "Reading file " << fname << "." << endl;
}

void FileIOUtils::NotifyFileWriting(char* fname)
{
  cout << "Writing file " << fname << "." << endl;
}

bool FileIOUtils::LineIsCommented(char* data) 
{
  bool move_on = false;

  if (data[0] == '#') move_on = true;       
  else if (data[0] == '\n') move_on = true; 

  return move_on;
}

char* FileIOUtils::StrUpp(char* string) 
{
  char* convert;

  convert = string;
  do { *convert = toupper((unsigned char)*convert); } while (*convert++);
  
  return string;
} 

char* FileIOUtils::StrLwr(char *string) 
{
  char *convert;

  convert = string;
  do { *convert = tolower((unsigned char)*convert); } while (*convert++);
  
  return string;

} 

