#include <iostream>
#include <stdio.h>

#include "simbox.h"
#include "exception.h"

using namespace std; 

int main(int argc, char** argv) {

  // setup simulation box
  Simbox* simbox;
  try {
    simbox = new Simbox();
  }
  catch (Exception* e) {
    e->Report();
    delete e;
    return 1;      
  }

  // run the simulation
  try {
    simbox->RunSimulation();
  }
  catch (Exception* e) {
    e->Report();
    delete e;
    return 1;
  }

  //clean up
  delete simbox;

  // no errors, return success
  return 0;
}
