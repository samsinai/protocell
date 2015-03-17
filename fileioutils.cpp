#include <iostream>
#include <fstream>
#include <fcntl.h>
<<<<<<< HEAD
#include <unistd.h>
=======
>>>>>>> 1e14aa9cd361d15fa8e0e0f7ebbd0b81c674e323
#include "exception.h"
#include "fileioutils.h"

using namespace std;

<<<<<<< HEAD
/* Find out if a file is binary or text stream */
// taken from PumMa
bool FileIOUtils::FileIsAscii(char* file)
{
  int i, ifile, in;
  char buf[512];
  
  ifile = open(file, 0);
  if(ifile < 0) 
    {
      throw new FileNotFoundException(file);

    }
  in = read(ifile, buf, 512);
  if(in == 0)
    {
      cout << "WARNING: " << file << " is an empty file" << endl;
      return -1;
    }
  for(i=0; i < in; i++) if(buf[i]&0200)
    {
      /* BINARY FILE */
      return false;
    }
  /* ASCII FILE */
  return true;
}

=======
>>>>>>> 1e14aa9cd361d15fa8e0e0f7ebbd0b81c674e323
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

