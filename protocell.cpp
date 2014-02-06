#include <fstream>
#include <iostream>
#include <stdio.h>

#include "protocell.h"
#include "mytypes.h"
#include "random.h"

Protocell::Protocell(long myid, long idparent) {
  protocellid = myid;
  parentid = idparent;
  Reset();
}

Protocell::~Protocell(void) {
  protocellcontents.clear();
}

void Protocell::Reset(void) {
  for (int i=0; i<NRMOLTYPES; i++) {
    nrs[i] = 0;
    initialnrs[i] = -1;
  }
  protocellcontents.clear();
  return;
}

void Protocell::AddMolecule(int n) {
  protocellcontents.push_back(n);
  nrs[n]++;
  return;
}

void Protocell::AddMoleculeFront(int n) {
  protocellcontents.push_front(n);
  nrs[n]++;
  return;
}

int Protocell::RemoveMoleculeAt(int n) {
  int q = protocellcontents[n];
  protocellcontents.erase(protocellcontents.begin()+n);
  nrs[q]--; 
  return q;
}

int Protocell::RemoveRandomMolecule(double rnd) {
  // rnd should be random number between 0 and 1
  int n = (int)(GetTotalNumberOfMolecules()*rnd);

  int q = RemoveMoleculeAt(n);
  return q;
}

int Protocell::RemoveMoleculeBack(void) {
  int q = protocellcontents.back();
  protocellcontents.pop_back();
  nrs[q]--; 
  return q;
}

Protocell* Protocell::SplitIntoTwo(double rnd, int& seed, int splitmethod) {
  // rnd should be random number between 0 and 1
  Protocell *newprotocell;
  int n = protocellcontents.size();

  if (n > 1) {
    newprotocell = new Protocell(-1,protocellid);
    int m;
    if (splitmethod == DIVIDE_LINEARPERTYPE) {
      int splitok = 0;
      while (splitok == 0) {
	for (int i=0; i<NRMOLTYPES; i++) {
	  int k = (int)(Random::Uniform(seed)*(nrs[i]+1));
	  int j = GetTotalNumberOfMolecules() - 1;
	  while ((k > 0) && (j > -1)) {
	    if (protocellcontents[j] == i) {
	      int q = RemoveMoleculeAt(j);
	      newprotocell->AddMoleculeFront(q);
	      k--;
	    }
	    j--;
	  }
	}
	if ((GetTotalNumberOfMolecules() > 0) && (newprotocell->GetTotalNumberOfMolecules() > 0)) {
	  splitok = 1;
	}
	else {
	  for (int j=0; j<newprotocell->GetTotalNumberOfMolecules(); j++) {
	      int q = newprotocell->RemoveMoleculeBack();
	      AddMoleculeFront(q);
	  }
	}
      }
    }
    else {
      if (splitmethod == DIVIDE_RANDOM) {
	vector<double> v(n,0);
	double d=1.0;
	m = n;
	for (int i=1; i<=(n/2); i++) {
	  d *= ((double)m)/i;
	  m--;
	  v[i]=d; 
	  v[n-i]=d;
	}
	for (int i=2; i<n; i++) {
	  v[i] += v[i-1];
	}
	d = v[n-1]*rnd;
	m = 1;
	while (d > v[m]) {
	  m++;
	}
      }
      else if (splitmethod == DIVIDE_EQUAL) {
	m = n/2;
      }
      else { //(splitmethod == DIVIDE_LINEAR)
	m = 1+(int)((n-1)*rnd);
      }
      
      for (int i=0; i<m; i++) {
	int q = RemoveMoleculeAt((int)(Random::Uniform(seed)*(n--)));
	newprotocell->AddMoleculeFront(q);
      }
    }
  }
  else {
    newprotocell = NULL;
  }

  return newprotocell;
}


int Protocell::GetTotalNumberOfMolecules() {
  int totnr = 0;
  for (int i=0; i<NRMOLTYPES; i++) {
    totnr += nrs[i];
  }
  return totnr;
}

int Protocell::GetNumberOfMoleculesOfType(int n) {
  return nrs[n];
}

int Protocell::GetInitialTotalNumberOfMolecules() {
  int totnr = 0;
  for (int i=0; i<NRMOLTYPES; i++) {
    totnr += initialnrs[i];
  }
  return totnr;
}

int Protocell::GetInitialNumberOfMoleculesOfType(int n) {
  return initialnrs[n];
}

void Protocell::SetInitialNumberOfMolecules(void) {
  for (int i=0; i<NRMOLTYPES; i++) {
    initialnrs[i] = nrs[i];
  }

  return; 
}

void Protocell::ReportComposition(void) {

  cout << "A: " << nrs[AMOL];
  cout << "\tB: " << nrs[BMOL];
  cout << "\tC: " << nrs[CMOL];
  cout << "\tD: " << nrs[DMOL];
  cout << "\tE: " << nrs[EMOL];
  cout << "\tF: " << nrs[FMOL];
  cout << endl;

  return;
}

