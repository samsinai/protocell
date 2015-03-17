#ifndef PROTOCELL
#define PROTOCELL

#include <stdio.h>
#include <deque>
#include <vector>
#include "mytypes.h"

using namespace std;

class Protocell {
 public:
  Protocell(long, long);
  ~Protocell();

  void Reset(void);
  void AddMolecule(int);
  void AddMoleculeFront(int);
  int RemoveMoleculeAt(int);
  int RemoveMoleculeBack(void);
  int RemoveRandomMolecule(double);
  Protocell* SplitIntoTwo(double, int&, int); 
  int GetTotalNumberOfMolecules(void);
  int GetNumberOfMoleculesOfType(int);
  int GetInitialTotalNumberOfMolecules(void);
  int GetInitialNumberOfMoleculesOfType(int);
  void SetInitialNumberOfMolecules(void);
  void ReportComposition(void);
  
  long protocellid;
  long parentid;

  deque<int> protocellcontents;

 private:
  int nrs[NRMOLTYPES];
  int initialnrs[NRMOLTYPES];
};

#endif
