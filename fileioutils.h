/*
** fielioutils.h
** 
** Made by Crijns
** 
** Started on  Wed Jun  4 17:03:54 2008 Crijns
** Last update Wed Jun  4 17:03:54 2008 Crijns
*/

#ifndef   	FILEIOUTILS_H_
# define   	FILEIOUTILS_H_


class FileIOUtils
{
 public:
<<<<<<< HEAD
    /*
   * Checks whether the file supplied is in ascii or binary format.
   */ 
  static bool FileIsAscii(char*);
=======
>>>>>>> 1e14aa9cd361d15fa8e0e0f7ebbd0b81c674e323
  /*
   * Checks whether the data supplied represents a comment, i.e.,
   * whether it starts with '#'.
   */
  static bool LineIsCommented(char*);
  /*
   * Prints a message stating the file that is currently read.
   */
  static void NotifyFileReading(char*);
  /*
   * Prints a message stating the file that is currently written.
   */
  static void NotifyFileWriting(char*);
  /*
   * Convert string to upper case
   */
  static char* StrUpp(char*);
  /*
   * Convert string to lower case
   */
  static char* StrLwr(char*);
};

#endif 	    /* !FIELIOUTILS_H_ */
