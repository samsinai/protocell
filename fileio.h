#ifndef FILEIO
#define FILEIO

using namespace std;

class Simbox;

/*
 * This class provides static methods that perform file reading and
 * writing tasks.
 */
class FileIO
{
 public:
  /*
   * Reads the input file and sets the DpdConfig object supplied acordingly.
   */
  static void ReadInputFile(Simbox*);
  static void DumpSystem(Simbox*, long);
  static void ReadSystem(Simbox*, long);

 private:
  static const int LINE_LENGTH = 200;
  static const int NAME_LENGTH = 80;

  /*
   * Reads an input argument with the name supplied from the open
   * file. The type of argument is passed by int - the value is stored
   * in the var to whic the void pointer points.
   */
  static void ReadInput(FILE*, string, int, void*);

  // for input file parsing - types
  static const int INPINT = 0;
  static const int INPREAL = 1;
  static const int INPCHAR = 2;
  static const int INPLONG = 3;

};
#endif
