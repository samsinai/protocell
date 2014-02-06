#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>

#include "fileio.h"
#include "simbox.h"
#include "exception.h"
#include "fileioutils.h"
#include "protocell.h"


void FileIO::ReadInputFile(Simbox* simbox)
{
  FILE* f;
  int nrlines;
  char line[LINE_LENGTH];

  // open parameter file
  char fname[NAME_LENGTH];
  sprintf(fname, "input.%s", "pci");
  if (simbox->verbose) {
    FileIOUtils::NotifyFileReading(fname);
  }
  f = fopen(fname, "r");
  if (f == NULL) {
    throw new FileNotFoundException(fname);
  }

  // find out nr of lines??
  nrlines = 0;
  while (fgets(line, LINE_LENGTH, f) != NULL) 
    { 
      if (FileIOUtils::LineIsCommented(line)) continue;
      ++nrlines;
    }   
  rewind(f);

  // integrator options
  ReadInput(f, "it", INPLONG, &(simbox->start_it));
  ReadInput(f, "max_it", INPLONG, &(simbox->n_iter));
  ReadInput(f, "seed", INPINT, &(simbox->seed));
  ReadInput(f, "log", INPINT, &(simbox->logreactions));

  // reaction rates
  ReadInput(f, "aa", INPREAL, &(simbox->prob_aA));
  ReadInput(f, "ac", INPREAL, &(simbox->prob_aC));
  ReadInput(f, "ad", INPREAL, &(simbox->prob_aD));
  ReadInput(f, "ae", INPREAL, &(simbox->prob_aE));
  ReadInput(f, "af", INPREAL, &(simbox->prob_aF));
  ReadInput(f, "k_s", INPREAL, &(simbox->prob_s));
  ReadInput(f, "q", INPREAL, &(simbox->prob_q));
  ReadInput(f, "k_phi", INPREAL, &(simbox->prob_phi));
  ReadInput(f, "k_theta", INPREAL, &(simbox->prob_theta));
  ReadInput(f, "k_z", INPREAL, &(simbox->prob_z));
  ReadInput(f, "z0", INPREAL, &(simbox->conc_z0));
  ReadInput(f, "k_d0", INPREAL, &(simbox->prob_d0));
  ReadInput(f, "k_d1", INPREAL, &(simbox->prob_d1));
  ReadInput(f, "m", INPINT, &(simbox->maxprotocellsize));
  ReadInput(f, "vol", INPREAL, &(simbox->volume));
  ReadInput(f, "revive", INPINT, &(simbox->revive));
  ReadInput(f, "modelnr", INPINT, &(simbox->modelnr));
  ReadInput(f, "splittype", INPINT, &(simbox->splittype));
  ReadInput(f, "intcomp", INPINT, &(simbox->intcomp));
  ReadInput(f, "divrep", INPREAL, &(simbox->divrep));

  // file dumping and reporting
  ReadInput(f, "erep", INPLONG, &(simbox->erep));
  ReadInput(f, "crep", INPLONG, &(simbox->crep));
  ReadInput(f, "timebrep", INPLONG, &(simbox->timebasedrep));

  // collecing data
  ReadInput(f, "avgrep", INPLONG, &(simbox->avgrep));
  ReadInput(f, "avgstart", INPLONG, &(simbox->avgstart));

  // initial composition
  ReadInput(f, "conc", INPREAL, &(simbox->conc));
  ReadInput(f, "maxnr", INPINT, &(simbox->maxnrprotocells));

  fclose(f);
}

void FileIO::ReadInput(FILE* f, string option, int type, void* optionVar)
{
  int linenr = 0;
  char line[LINE_LENGTH];
  char arg[NAME_LENGTH], val[NAME_LENGTH];
  bool optionEncounter = false;
  char *cptr;
  int *iptr;
  long *lptr;
  double *rptr;

  while (fgets(line, LINE_LENGTH, f) != NULL) 
    { 
      if (FileIOUtils::LineIsCommented(line)) { linenr ++;  continue; }
      sscanf(line, "%s %s", arg, val); 

      FileIOUtils::StrLwr(arg);
      FileIOUtils::StrLwr(val);
      if (!option.compare(arg) && !optionEncounter)
        {
	  //cout  << "reading " << arg << endl;

          optionEncounter = true;
          switch (type)
            {
            case (INPINT):
              iptr = (int*)optionVar; 
              *iptr = atoi(val);
              break;
            case (INPLONG):
              lptr = (long*)optionVar; 
              *lptr = atol(val);
              break;
            case (INPREAL):
              rptr = (double*)optionVar; 
              *rptr = atof(val);
              break;
            case (INPCHAR):
              cptr = (char*)optionVar;
              strcpy(cptr,val);
              break;
            default:
              throw new Exception("Unknown input type specified");
            }
        }
      else if (!option.compare(arg) && optionEncounter)
        {
          throw new 
            Exception("Encountered option %s for a second time on line %d", 
                      option.c_str(), linenr);
        }
      linenr ++;
    }
  rewind(f);
}

  
void FileIO::DumpSystem(Simbox* simbox, long it) {
  fstream of;
  
  char fname[NAME_LENGTH];
  sprintf(fname, "IT%ld.pcs", it);
  of.open(fname, ios::out);
  if(!of.is_open()) {
    throw new Exception("Could not open file '%s'.", fname);
  }
  of << "<system time=\"" << simbox->time;
  of << "\">" << endl;
  of << "<z0> ";
  of << simbox->nr_zs;
  of << " </z0>" << endl;
  for (int mi=0; mi<simbox->nrprotocells; mi++) {
    Protocell *m;
    m = simbox->protocells[simbox->usednrs[mi]];
    of << "  <protocell ";
    of << "id=\"" << m->protocellid;
    of << "\" parentid=\"" << m->parentid;
    of << "\"> "; // << endl << "     ";
    //int totj = 0;
    for (int j=0; j<m->GetTotalNumberOfMolecules(); j++) {
      int i = m->protocellcontents[j];
      //totj++;
      if ((i >= 0) && (i < NRMOLTYPES)) {
	of << moltypenames[i];
      }
      else {
	of << "?";
      }
      //if ((totj%10) == 0) {
      //  of << " ";
      //}
      //if (((totj%60) == 0) && (totj!=m->GetTotalNumberOfMolecules())) {
      //  of << endl << "    ";
      //}
      //}
    }
    //of << endl;
    of << " </protocell>" << endl;
  }

  of << "</system>" << endl;

  of.close();

  return;
}


void FileIO::ReadSystem(Simbox* simbox, long it) {
  fstream inf;
  
  char fname[NAME_LENGTH];
  sprintf(fname, "IT%ld.pcs", it);
  inf.open(fname, ios::in);
  if(!inf.is_open()) {
    throw new Exception("Could not open file '%s'.", fname);
  }
  cout << "Reading configuration from file \"" << fname;
  cout << "\"" << endl;

  string tc;
  inf >> tc;
  if (tc.compare("<system")) {
    cout << "Found: " << tc << endl;
    throw new Exception("Expected file starting with: \"<system\"");
  }
  inf >> tc;
  if (tc.compare(0,6,"time=\"")) {
    cout << "Found: " << tc << endl;
    throw new Exception("Expected time=\"...\"");
  }
  char buf[20]="                   ";
  tc.copy(buf,tc.size()-8, 6);
  //cout << "qqq" << buf << "qqq" << endl; 
  simbox->time = atof(buf);
  simbox->nrprotocells = 0;
  bool finished = false;
  while (!inf.eof() && !finished) {
    inf >> tc;
    if (!tc.compare("</system>")) {
      finished = true;
    }
    else if (!tc.compare("<z0>")) {
      inf >> tc;
      simbox->nr_zs=atoi(tc.c_str());
      inf >> tc;
      if (tc.compare("</z0>")) {
	cout << "Found: " << tc << endl;
	throw new Exception("Expected: \"</z0>\"");
      }
    }
    else if (!tc.compare("<protocell")) {
      //cout << "Reading a protocell" << endl;
      inf >> tc;
      Protocell *newprotocell;
      newprotocell = new Protocell(-1,-1);
      if (!tc.compare(0,4,"id=\"")) {
	char buf[20]="                   ";
	tc.copy(buf,tc.size()-5, 4);
	//cout << "ppp" << buf << "pp" << endl;
	newprotocell->protocellid = atol(buf);
	if ((newprotocell->protocellid) >= simbox->nextuid) {
	  simbox->nextuid = newprotocell->protocellid + 1;
	}
      }
      else {
	cout << tc << endl;
	throw new Exception("Protocell id expected");	
      }
      inf >> tc;
      if (!tc.compare(0,10,"parentid=\"")) {
	char buf[20]="                   ";
	tc.copy(buf,tc.size()-12, 10);
	//cout << "rrr" << buf << "rr" << endl;
	newprotocell->parentid = atol(buf);
      }
      else {
	cout << tc << endl;
	throw new Exception("Protocell parentid expected");	
      }
      inf >> tc;
      while (!inf.eof() && (tc.compare("</protocell>"))) {
	for(unsigned int i=0; i<tc.length(); i++) {
	  switch (tc.at(i)) {
	  case (AMOLNAME):
	    newprotocell->AddMolecule(AMOL);
	    break;
	  case (BMOLNAME):
	    newprotocell->AddMolecule(BMOL);
	    break;
	  case (CMOLNAME):
	    newprotocell->AddMolecule(CMOL);
	    break;
	  case (DMOLNAME):
	    newprotocell->AddMolecule(DMOL);
	    break;
 	  case (EMOLNAME):
 	    newprotocell->AddMolecule(EMOL);
 	    break;
 	  case (FMOLNAME):
 	    newprotocell->AddMolecule(FMOL);
 	    break;
	  default:
	    cout << tc.at(i) << endl;
	    throw new Exception("Unknown molecule type");
	  }
	}
	inf >> tc;
      }  
      simbox->AddProtocell(newprotocell,0);
    }
    else {
      cout  << endl << tc << endl;
      throw new Exception("Unexpected keyword");
    }
  }
  if (inf.eof()) {
    cout << endl;
    throw new Exception("Unexpected end of file\n");    
  }
  
  return;
}
