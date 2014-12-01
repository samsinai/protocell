#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>

#include "exception.h"
#include "simbox.h"
#include "fileio.h"
#include "random.h"
#include "protocell.h"
#include "mytypes.h"


Simbox::Simbox(void) {
  time = 0.0;
  start_it = -1;
  nextuid = 0;
  timebasedrep = 0;
  erep = 100;
  crep = -1;
  nrprotocells = -1;
  maxnrprotocells = 100000;
  n_iter = 100;
  conc = -1.0;
  summedpropensities = 0.0;
  logreactions = 0;
  avgrep = -1;
  avgstart = -1;
  splittype = DIVIDE_RANDOM;
  intcomp = 0;
  divrep = 1.0;
  verbose = 1;

  prob_aA = -1.0;
  prob_aC = -1.0;
  prob_aD = -1.0;
  prob_aE = -1.0;
  prob_aF = -1.0;
  prob_s = -1.0;
  prob_q = -1.0;
  prob_phi = -1.0;
  prob_theta = 0.0;
  prob_z = -1.0;
  conc_z0 = -1.0;
  prob_d0 = 0.0;
  prob_d1 = 0.0;
  maxprotocellsize= -1;
  volume = -1;
  nr_zs = 0;

  revive = -1;
  nrBcells = 0;
  nrrevivals = 0;
  modelnr = -1;

  FileIO::ReadInputFile(this);

  if (start_it == -1) {throw new Exception("start_it not (or invalidly) specified"); }
  if (prob_aA < 0.0) { throw new Exception("aA not (or invalidly) specified"); }
  if (prob_aC < 0.0) { throw new Exception("aC not (or invalidly) specified"); }
  if (prob_aD < 0.0) { throw new Exception("aD not (or invalidly) specified"); }
  if (prob_aE < 0.0) { throw new Exception("aE not (or invalidly) specified"); }
  if (prob_aF < 0.0) { throw new Exception("aF not (or invalidly) specified"); }
  if (prob_s < 0.0) { throw new Exception("k_s not (or invalidly) specified"); }
  if (prob_q < 0.0) { throw new Exception("q not (or invalidly) specified"); }
  if (prob_z < 0.0) { throw new Exception("k_z not (or invalidly) specified"); }
  if (conc_z0 < 0.0) { throw new Exception("z0 not (or invalidly) specified"); }
  if (prob_phi < 0.0) { throw new Exception("k_phi not (or invalidly) specified"); }
  if (prob_theta < 0.0) { throw new Exception("k_theta invalidly specified"); }
  if (maxprotocellsize < 0) { throw new Exception("m not (or invalidly) specified"); }
  if (volume < 0.0) { throw new Exception("vol not (or invalidly) specified"); }
  if (modelnr < 1) { throw new Exception("modelnr not (or invalidly) specified"); }

  invvol = 1.0/volume;

  CalculateMesoscopicRateConstants();

  usednrs.reserve((size_t)(maxnrprotocells));
  availablenrs.resize(maxnrprotocells);
  for (int i=0; i<maxnrprotocells; i++) {
    availablenrs[i]=maxnrprotocells-1-i;
  }
  snrinfo.resize(maxnrprotocells, -1);
  nrpropsperprotocell = 9;
  propensities.resize(nrpropsperprotocell*maxnrprotocells, 0.0);

  protocells.resize(maxnrprotocells);
  for (int i=0; i<maxnrprotocells; i++) {
    protocells[i] = NULL;
  }
  for (int i=0; i<NRMOLTYPES; i++) {
    totnrmolspertype[i] = 0;
  }

  cout << "************************************" << endl;
  cout << "*          ProtoCell v0.2          *" << endl;
  cout << "*          A.J. Markvoort          *" << endl;
  cout << "************************************" << endl;
  cout << endl;

  cout << "  prob_aA:     " << prob_aA << endl;
  cout << "  prob_aC:     " << prob_aC << endl;
  cout << "  prob_aD:     " << prob_aD << endl;
  cout << "  prob_aE:     " << prob_aE << endl;
  cout << "  prob_aF:     " << prob_aF << endl;
  cout << "  prob_s:     " << prob_s << endl;
  cout << "  prob_q:     " << prob_q << endl;
  cout << "  prob_phi:   " << prob_phi << endl;
  cout << "  prob_theta: " << prob_theta << endl;
  cout << "  prob_d0:    " << prob_d0 << endl;
  cout << "  prob_d1:    " << prob_d1 << endl;
  cout << "  prob_z:     " << prob_z << endl;
  cout << "  conc_z0 :   " << conc_z0 << endl;
  cout << "  m:          " << maxprotocellsize << endl;

  cout << endl;
  cout << "Number of iterations       " << n_iter << endl;
  cout << "Random seed                " << seed << endl;
  if (timebasedrep == 0) {
    cout << "Report frequency (erep)    " << erep << endl;
  }
  else {
    nexttimebasedrep = 0.0;
    cout << "Time based report" << endl;
  }
  if (crep < 0) {
    cout << "Config. dump off" << endl;
  }
  else {
    cout << "Config. dump freq. (crep)  " << crep << endl;
  }

  cout << "Split type: ";
  if (splittype == DIVIDE_RANDOM) {
    cout << "DIVIDE_RANDOM";
  }
  else if (splittype == DIVIDE_EQUAL) {
    cout << "DIVIDE_EQUAL";
  }
  else if (splittype == DIVIDE_LINEAR) {
    cout << "DIVIDE_LINEAR";
  }
  else if (splittype == DIVIDE_LINEARPERTYPE) {
    cout << "DIVIDE_LINEARPERTYPE";
  }
  else if (splittype == DIVIDE_EQUALPERTYPE) {
    cout << "DIVIDE_EQUALPERTYPE";
  }
  else {
    cout << "DIVIDE_MANY";
  }
  cout << endl;
  cout << "divrep: " << divrep << endl;

  if (intcomp == 1) {
    cout << "Internal competition\n" << endl;
  }
  else {
    cout << "No internal competition\n" << endl;
  }



  cout << "Starting from iteration " << start_it << " at time point " << time << endl;

  FileIO::ReadSystem(this, start_it);
  
  cout << "Initial number of protocells = " << nrprotocells << endl;
  //cout << "Concentration:             " << conc << " mol/l" << endl;
  cout << "Volume:                    " << volume << " l" << endl;
  cout << endl;

  conc_z0 *= 6.0221415e23*volume;

  if (avgrep == -1) {
    cout << "Not collecting any average data." << endl;
    avgrepnext = n_iter+10;
  }
  else {
    if (avgstart == -1) {
      avgstart = start_it;
    }
    avgrepnext = avgstart;
    cout << "Data collect freq. (avgrep)" << avgrep << endl;
    cout << "Data coll. start (avgstart)" << avgstart << endl;
  }
  cout << endl;

  if (revive > -1) {
    if (revive < NRMOLTYPES) {
      cout << "Will create new cell with single nolecule of type " << revive << "when last molecule of that type disappeared" << endl;
      cout << endl;
    }
    else {
      throw new Exception("Too large number for revive (%d). Should be smaller than number of molecule types", revive);
    }
  }

  for (int i=0; i<NRMOLTYPES; i++) {
    nrcreations[i] = 0;
    nrremovals[i] = 0;
  }
  nrprotocelldivisions = 0;

  return;
}


Simbox::~Simbox(void) {
  for (int i=0; i<nrprotocells; i++) {
    delete protocells[usednrs[i]];
  }
  protocells.clear();
  availablenrs.clear();
  usednrs.clear();
  snrinfo.clear();
  propensities.clear();
}


void Simbox::CalculateMesoscopicRateConstants(void) {
  // conversion to mesoscopic rate constants
  prob_s *= invvol/6.0221415e23; 
  //prob_z0 *= 6.0221415e23*volume;

  return;
}


void Simbox::AddProtocell(Protocell *newprotocell, int resultofdivision) {

  int newprotocellnr = availablenrs.back();
  availablenrs.pop_back();
  for (int i=0; i<nrpropsperprotocell; i++) {
    propensities[nrpropsperprotocell*newprotocellnr+i] = 0;
  }
  protocells[newprotocellnr] = newprotocell;
  summedpropensities += calculate_propensities_protocell(newprotocellnr);

  if (resultofdivision == 0) {
    for (int i=0; i<NRMOLTYPES; i++) {
      totnrmolspertype[i] += newprotocell->GetNumberOfMoleculesOfType(i);
    }
  }

  if ((newprotocell->GetTotalNumberOfMolecules() - newprotocell->GetNumberOfMoleculesOfType(BMOL)) == 0) { nrBcells++; }

  snrinfo[newprotocellnr]=usednrs.size();
  usednrs.push_back(newprotocellnr);
  nrprotocells++;

  return;
}


void Simbox::RemoveProtocell(int nr) {
  Protocell *m;
  //cout << "rrr " << nr << " " << usednrs.size() << endl;

  m = protocells[nr];

  //cout << "111" << endl;
  if ((m->GetTotalNumberOfMolecules() - m->GetNumberOfMoleculesOfType(BMOL)) == 0) { nrBcells--; }
  //cout << "222" << endl;

  for (int i=0; i<NRMOLTYPES; i++) {
    totnrmolspertype[i] -= m->GetNumberOfMoleculesOfType(i);
  }
  //cout << "333" << endl;
  for (int i=0; i<nrpropsperprotocell; i++) {
    summedpropensities -= propensities[nrpropsperprotocell*nr+i];
    propensities[nrpropsperprotocell*nr+i] = 0.0;
  }
  //cout << "444" << endl;
  delete m;
  protocells[nr] = NULL;
  nrprotocells--;

  //cout << "555" << endl;
  snrinfo[usednrs.back()]=snrinfo[nr];
  usednrs[snrinfo[nr]]=usednrs.back();
  snrinfo[nr]=-1;
  usednrs.pop_back();
  availablenrs.push_back(nr);
  //cout << "666" << endl;

  return;
}


void Simbox::DivideProtocell(double rnd, int cellnr, ostream& s, int it) {
  Protocell *m;

  m = protocells[cellnr];
  if ((m->GetTotalNumberOfMolecules() - m->GetNumberOfMoleculesOfType(BMOL)) == 0) { nrBcells--; }
  Protocell *mnew;
  //cout << "splitting at time " << time << endl;
  if (m->GetTotalNumberOfMolecules() > 1) {
    double myinitialfrac = -1.0;
    if (it > avgstart) {

      int mynrmols = m->GetTotalNumberOfMolecules();
      int nrnonB = mynrmols - m->GetNumberOfMoleculesOfType(1);
      double myfract = (100.0*nrnonB)/mynrmols;
      avgcol_fractionbeforedivision[(int)myfract] += 1;
      
      //if (nrnonB == 0) { 
      //	avgcol_nrcellsbeforedivision0C += 1;
      //	avgcol_cellsizebeforedivision0C += mynrmols;
      //}
      //else if (nrnonB == 1) { 
      //	avgcol_nrcellsbeforedivision1C += 1;
      //	avgcol_cellsizebeforedivision1C += mynrmols;
      //}
      //else { 
      //	avgcol_nrcellsbeforedivisionmC += 1;
      //	avgcol_cellsizebeforedivisionmC += mynrmols;
      //}
      //if (nrnonB == mynrmols) { 
      //	avgcol_nrcellsbeforedivisionpC += 1;
      //	avgcol_cellsizebeforedivisionpC += mynrmols;
      //}
      
      int myinitialnrmols = m->GetInitialTotalNumberOfMolecules();
      if (myinitialnrmols > 0) {
	int initialnrnonB = myinitialnrmols - m->GetInitialNumberOfMoleculesOfType(1);
	myinitialfrac = (100.0*initialnrnonB)/myinitialnrmols;
      }
    }

    if (splittype == DIVIDE_MANY) {
      while (m->GetTotalNumberOfMolecules() > 1) {
	mnew = new Protocell(nextuid++, m->protocellid);
	int q = m->RemoveMoleculeBack();
	mnew->AddMoleculeFront(q);
	AddProtocell(mnew,1);
	mnew->SetInitialNumberOfMolecules();
      }
      if (logreactions) {
	s << time << " split protocell " << m->protocellid;
	s << " creating many";
	s << endl;
      }
    }
    else {
      mnew = m->SplitIntoTwo(rnd, seed, splittype);
      mnew->protocellid = nextuid++;

      if (it > avgstart) {
	int mynrmols = m->GetTotalNumberOfMolecules();
	int nrnonB = mynrmols - m->GetNumberOfMoleculesOfType(1);
	double myfract = (100.0*nrnonB)/mynrmols;
	avgcol_fractionafterdivision[(int)myfract] += 1;
	
	//if (nrnonB == 0) { 
	//  avgcol_nrcellsafterdivision0C += 1;
	//  avgcol_cellsizeafterdivision0C += mynrmols;
	//}
	//else if (nrnonB == 1) { 
	//  avgcol_nrcellsafterdivision1C += 1;
	//  avgcol_cellsizeafterdivision1C += mynrmols;
	//}
	//else { 
	//  avgcol_nrcellsafterdivisionmC += 1;
	//  avgcol_cellsizeafterdivisionmC += mynrmols;
	//}
	//if (nrnonB == mynrmols) { 
	//  avgcol_nrcellsafterdivisionpC += 1;
	//  avgcol_cellsizeafterdivisionpC += mynrmols;
	//}
	
	mynrmols = mnew->GetTotalNumberOfMolecules();
	nrnonB = mynrmols - mnew->GetNumberOfMoleculesOfType(1);
	double myfract2 = (100.0*nrnonB)/mynrmols;
	avgcol_fractionafterdivision[(int)myfract2] += 1;
	
	//if (nrnonB == 0) { 
	//  avgcol_nrcellsafterdivision0C += 1;
	//  avgcol_cellsizeafterdivision0C += mynrmols;
	//}
	//else if (nrnonB == 1) { 
	//  avgcol_nrcellsafterdivision1C += 1;
	//  avgcol_cellsizeafterdivision1C += mynrmols;
	//}
	//else { 
	//  avgcol_nrcellsafterdivisionmC += 1;
	//  avgcol_cellsizeafterdivisionmC += mynrmols;
	//}
	//if (nrnonB == mynrmols) { 
	//  avgcol_nrcellsafterdivisionpC += 1;
	//  avgcol_cellsizeafterdivisionpC += mynrmols;
	//}
	
	if (myinitialfrac > -0.5) {
	  avgcol_fractioninitvsfinal[(int)myinitialfrac][(int)myfract] += 1;
	  avgcol_fractioninitvsfinal[(int)myinitialfrac][(int)myfract2] += 1;
	}
      }

      AddProtocell(mnew,1);
      mnew->SetInitialNumberOfMolecules();

      if (logreactions) {
	s << time << " split protocell " << m->protocellid;
	s << " creating " << mnew->protocellid << " ";
	for (int i=0; i< mnew->GetTotalNumberOfMolecules(); i++) {
	  s << moltypenames[mnew->protocellcontents[i]];
	}
	s << endl;
      }
    }
    summedpropensities -= get_storedtotalpropensity_protocell(cellnr);
    summedpropensities += calculate_propensities_protocell(cellnr);
    m->SetInitialNumberOfMolecules();
    nrprotocelldivisions++;
    if ((m->GetTotalNumberOfMolecules() - m->GetNumberOfMoleculesOfType(BMOL)) == 0) { nrBcells++; }

  //cout << "done splitting" << endl;
  }  
  return;
}


int Simbox::perform_reaction(int nr, double rnd, ostream& s, int it) {
  int cellnr = nr/nrpropsperprotocell;
  int reacnr = nr%nrpropsperprotocell;
  Protocell *m;
  //cout << "ccc " << cellnr << " " << snrinfo[cellnr] << " " << rnd <<endl; 
  m = protocells[cellnr];
  //cout << "Performing reaction " << reacnr << " on protocell " << m->protocellid << " (" << cellnr << ")"<< endl;

  if (reacnr < 6) {
    m->AddMolecule(reacnr);
    nr_zs--;
    summedpropensities -= get_storedtotalpropensity_protocell(cellnr);
    summedpropensities += calculate_propensities_protocell(cellnr);
    totnrmolspertype[reacnr]++;
    nrcreations[reacnr]++;
    if (logreactions) {
      s << time << " add " << moltypenames[reacnr] << " to protocell " << m->protocellid << endl;
    }
  }
  else if (reacnr == 6) {
    // protocell dies
    if (logreactions) {
      s << time << " remove protocell " << m->protocellid << endl;
    }
    RemoveProtocell(cellnr);
    m = NULL;
    //cout << "ttt" << endl;
  }
  else if (reacnr == 7) {
    // single molecule in protocell dies
    if (m->GetTotalNumberOfMolecules() == 1) {
      int q = m->protocellcontents.back();
      nrremovals[q]++;
      RemoveProtocell(cellnr);
      m = NULL;
      if (logreactions) {
	s << time << " remove last molecule from protocell " << m->protocellid << endl;
      }
    }
    else {
      int q = m->RemoveRandomMolecule(rnd);
      summedpropensities -= get_storedtotalpropensity_protocell(cellnr);
      summedpropensities += calculate_propensities_protocell(cellnr);
      totnrmolspertype[q]--;
      nrremovals[q]++;
      if (logreactions) {
	s << time << " remove " << moltypenames[q] << " from protocell " << m->protocellid << endl;
      }
    }
  }
  else { // reacnr == 8
    //protocell divides 
    DivideProtocell(rnd, cellnr, s, it);
  }
  //cout << "uuu" << endl;

  // division of protocell when maximum protocell size is reached 
  if ((m != NULL) && (m->GetTotalNumberOfMolecules() == maxprotocellsize)) {
    if (rnd < divrep) { // divide or remove one molecule from cell
      rnd = rnd/divrep;
      DivideProtocell(rnd, cellnr, s, it);
    }
    else {
      rnd = (rnd-divrep)/(1.0-divrep);
      int q = m->RemoveMoleculeAt((int)(rnd*m->GetTotalNumberOfMolecules()));
      summedpropensities -= get_storedtotalpropensity_protocell(cellnr);
      summedpropensities += calculate_propensities_protocell(cellnr);
      totnrmolspertype[q]--;
      nrremovals[q]++;
      if (logreactions) {
	s << time << " remove " << moltypenames[q] << " from protocell " << m->protocellid << " at division threshold" << endl;
      }
    }
  }
  //cout << "vvv" << endl;

  return reacnr;
}


double Simbox::get_storedtotalpropensity_protocell(int nr) {
  double deltaptot = 0.0;

  for (int i=0; i<nrpropsperprotocell; i++) {
    deltaptot += propensities[nrpropsperprotocell*nr+i];
  }

  return deltaptot;
}


double Simbox::calculate_propensities_protocell(int nr) {
  double deltaptot = 0.0;
  Protocell *m;
  m = protocells[nr];
  
  int nA = m->GetNumberOfMoleculesOfType(AMOL);
  int nB = m->GetNumberOfMoleculesOfType(BMOL);
  int nC = m->GetNumberOfMoleculesOfType(CMOL);
  int nD = m->GetNumberOfMoleculesOfType(DMOL);
  int nE = m->GetNumberOfMoleculesOfType(EMOL);
  int nF = m->GetNumberOfMoleculesOfType(FMOL);
  int ntot = nA+nB+nC+nD+nE+nF;

  int nAm1 = nA-1;
  if (nAm1 < 0) {
    nAm1 = 0;
  }

  if (modelnr == 1) { // model R2alpha
    double p0 = nr_zs*(1+nC*prob_aC)*prob_s;

    if (intcomp == 1) {
      p0 /= ntot;
    }
    propensities[nrpropsperprotocell*nr+0] = p0*(nA*(1+nAm1*prob_aA)*prob_q);
    propensities[nrpropsperprotocell*nr+1] = p0*(nA*(1+nAm1*prob_aA)*(1-prob_q) + 
						 nB*(1+nA*prob_aA) + 
						 nC*(1+nA*prob_aA)*(1-prob_q) +
						 nD*(1+nA*prob_aA)*(1-prob_q) +
						 nE*(1+nA*prob_aA)*(1-prob_q) +
						 nF*(1+nA*prob_aA)*(1-prob_q));
    propensities[nrpropsperprotocell*nr+2] = p0*(nC*(1+nA*prob_aA)*prob_q);
    propensities[nrpropsperprotocell*nr+3] = p0*(nD*(1+nA*prob_aA)*prob_q);
    propensities[nrpropsperprotocell*nr+4] = p0*(nE*(1+nA*prob_aA)*prob_q);
    propensities[nrpropsperprotocell*nr+5] = p0*(nF*(1+nA*prob_aA)*prob_q);
    propensities[nrpropsperprotocell*nr+6] = (1+prob_aE*nE)*prob_phi;
    propensities[nrpropsperprotocell*nr+7] = (1+prob_aF*nF)*ntot*prob_theta;
    propensities[nrpropsperprotocell*nr+8] = (1+prob_aD*nD)*(prob_d0+ntot*prob_d1);
    if (propensities[nrpropsperprotocell*nr+8] < 0.0) {
      propensities[nrpropsperprotocell*nr+8] = 0.0;
    }
  }
  else if (modelnr == 2) { // model R1 (C) and R2 (A)
    double mfAs = (1+(prob_aA-1)*(nA>1));  // selfamplification only when more than 1 A present
    double mfAo = (1+(prob_aA-1)*(nA>0));  // for amplification of others 1 A suffices
    double mfC = (1+(prob_aC-1)*(nC>0));
    double mfD = (1+(prob_aD-1)*(nD>0));
    double mfE = (1+(prob_aE-1)*(nE>0));
    double mfF = (1+(prob_aF-1)*(nF>0));
    double p0 = nr_zs*mfC*prob_s;
    if (intcomp == 1) {
      p0 /= ntot;
    }
    propensities[nrpropsperprotocell*nr+0] = p0*nA*mfAs*prob_q;
    propensities[nrpropsperprotocell*nr+1] = p0*(nA*mfAs*(1-prob_q) + 
						 nB*mfAo +
						 nC*mfAo*(1-prob_q) +
						 nD*mfAo*(1-prob_q) +
						 nE*mfAo*(1-prob_q) +
						 nF*mfAo*(1-prob_q));
    propensities[nrpropsperprotocell*nr+2] = p0*nC*mfAo*prob_q;
    propensities[nrpropsperprotocell*nr+3] = p0*nD*mfAo*prob_q;
    propensities[nrpropsperprotocell*nr+4] = p0*nE*mfAo*prob_q;
    propensities[nrpropsperprotocell*nr+5] = p0*nF*mfAo*prob_q;
    propensities[nrpropsperprotocell*nr+6] = mfE*prob_phi;
    propensities[nrpropsperprotocell*nr+7] = mfF*ntot*prob_theta;
    propensities[nrpropsperprotocell*nr+8] = mfD*(prob_d0+ntot*prob_d1);
    if (propensities[nrpropsperprotocell*nr+8] < 0.0) {
      propensities[nrpropsperprotocell*nr+8] = 0.0;
    }
  }

  for (int i=0; i<nrpropsperprotocell; i++) {
    deltaptot += propensities[nrpropsperprotocell*nr+i];
  }

  return deltaptot;
}


void Simbox::RunSimulation(void) {
  fstream of, lf;
  char fname[80];
  sprintf(fname, "output.%s", "pco");
  of.open(fname, ios::out);
  if(!of.is_open()) {
    throw new Exception("Could not open file %s.", fname);
  }

  if (logreactions) {
    sprintf(fname, "reactionlog.pcr");
    lf.open(fname, ios::out);
    if(!lf.is_open()) {
      throw new Exception("Could not open file %s.", fname);
    }
  }

  if (avgrep != -1) {
    InitAveragesCollection();
  }

  long crepnext;
  if (crep < 0) {
    crepnext = n_iter;
  }
  else {
    crepnext = start_it+crep; 
  }

  summedpropensities += prob_z*fabs(conc_z0-nr_zs);
  cout << "Total propensity at initial time point: " << summedpropensities << endl;
  ReportStatus(cout, start_it);
  ReportStatus(of, start_it);

  int nrselected = -1;
  long maxit = n_iter;
  long it = start_it+1; 
  while ((it <= maxit) && ((nrprotocells > 0) || (revive > -1))) {
    bool incit = true;

    if (revive > -1) {
      if (totnrmolspertype[revive] == 0) {
	Protocell *newprotocell;
	newprotocell = new Protocell(-1,-1);
	newprotocell->protocellid = nextuid++;
	newprotocell->AddMolecule(revive);
	AddProtocell(newprotocell,0);
	nrrevivals++;
      }
    }

    double prob_z0 = prob_z*fabs(conc_z0-nr_zs);

    // recalculate all propensities because they all depend on z0
    summedpropensities = 0.0;
    for (int i=0; i<nrprotocells; i++) {
      summedpropensities += calculate_propensities_protocell(usednrs[i]);
    }
    summedpropensities += prob_z0;

    double r = Random::Uniform(seed);
    r *= summedpropensities;
    
    // increase time
    double y = Random::Uniform(seed);
    double deltatime = -log(y)/summedpropensities;
    time += deltatime;
    
    //if (it>4990000) cout << "aaa " << it << endl;
    // Perform the selected reaction
    if (r < prob_z0) {
      if (conc_z0 > nr_zs) {
	nr_zs++;
      }
      else {
	nr_zs--;
      }
      incit = false;
      //if (it>4990000) cout << "bbb " << endl;
    }
    else {
      r -= prob_z0;

      nrselected = 0;
      while (r >= propensities[nrselected]) { // really inefficient now because it possibly also loops over loops of zeros 
	r -= propensities[nrselected];
	nrselected++;
      }
      //if (it>4990000) cout << "bbb " << nrselected << endl;
      perform_reaction(nrselected, r/propensities[nrselected], lf, it);
      //if (it>4990000) cout << "www" << endl;
    }
    //if (it>4990000) cout << "ccc " << endl;

    if (incit == true) {

      if (timebasedrep == 0) {
	if ((it%erep) == 0) {
	  ReportStatus(cout, it);
	  ReportStatus(of, it);
	}
      }
      else {
	if (time > nexttimebasedrep) {
	  ReportStatus(cout, it);
	  ReportStatus(of, it);
	  nexttimebasedrep = 1.01*time;
	}
      }
      //if (it>4990000) cout << "ddd " << endl;
      
      if (it >= avgrepnext) {
	UpdateAveragesCollection(it, deltatime);
	avgrepnext += avgrep;
      }
      //if (it>4990000) cout << "eee " << endl;
      
      if ((crep > -1) && (it >= crepnext)) {
	FileIO::DumpSystem(this, it);
	crepnext += crep;      
      }
      //if (it>4990000) cout << "fff " << endl;
      
      if (((int) nrprotocells) >= maxnrprotocells) {
	throw new Exception("Maximal number of protocells reached.");
      }
      //if (it>4990000) cout << "ggg " << endl;
      
      it++;
    }
  }
  it--;
  ReportStatus(cout, it);
  ReportStatus(of, it);
  if (nrprotocells == 0) {
    cout << "No protocells left at iteration " << it << endl;
  }

  of.close(); // close output file
  if (logreactions) {
    lf.close(); // close reaction log file
  }
  if (avgrep != -1) {
    ReportAveragesCollection();
  }

  FileIO::DumpSystem(this, it);

  return;
}


void Simbox::ReportStatus(ostream& s, long& it) {

  s << "IT " << it  << " " << time << " " << nrprotocells;
  s << " " << totnrmolspertype[AMOL] << " " << totnrmolspertype[BMOL];
  s << " " << totnrmolspertype[CMOL] << " " << totnrmolspertype[DMOL];
  s << " " << totnrmolspertype[EMOL] << " " << totnrmolspertype[FMOL];
  s << " " << nr_zs;
  s << " " << nrrevivals;

  int nCD = 0;
  int npureC = 0;
  int npureD = 0;
  int npureCD = 0;
  for (int i=0; i<nrprotocells; i++) {
    int mtot = protocells[usednrs[i]]->GetTotalNumberOfMolecules();
    int mC = protocells[usednrs[i]]->GetNumberOfMoleculesOfType(CMOL);
    int mD = protocells[usednrs[i]]->GetNumberOfMoleculesOfType(DMOL);
    if ((mC*mD) > 0) { nCD++; }
    if (mC == mtot) { npureC++; }
    else if (mD == mtot) { npureD++; }
    else if ((mC+mD) == mtot){ npureCD++; }
  }
  s << " " << nCD << " " << npureC << " " << npureD << " " << npureCD;

  s << endl;

  return;
}


void Simbox::Check_propensity_consistency(void) {
  
  bool ok = true;
  for (int i=0; i<nrprotocells; i++) {
    if (protocells[usednrs[i]] == NULL) {
      ok = false;
    }
  }
  if (ok == false) {
    throw new Exception("Oops, protocell lost");
  }

  ok = true;
  for (int i=0; i<nrprotocells; i++) {
    if (snrinfo[usednrs[i]] == -1) {
      ok = false;
    }
  }
  if (ok == false) {
    throw new Exception("Oops, protocell lost");
  }

  //double checksummedpropensities = 0.0;
  //for (int i=0; i<nrprotocells; i++) {
  //  checksummedpropensities += calculate_propensities_protocell(usednrs[i]);
  //}
  //checksummedpropensities += prob_z0;
  //if (fabs(checksummedpropensities-summedpropensities) > 1e-3) {
  //  cout << "Checking total propensity:" << checksummedpropensities << " <> " << summedpropensities << endl;
  //  throw new Exception("Oops");
  //}

  return;
}


void Simbox::InitAveragesCollection(void) {
  nrhistupdates = 0;
 // nrABperCell(100, vector<int>(100));
// added by sam
  for (int i=0; i<101; i++){
    avgcol_fractionpresent[i] = 0.0;
    avgcol_fractionbeforedivision[i] = 0.0;
    avgcol_fractionafterdivision[i] = 0.0;

    for (int j=0; j<101; j++){
      avgcol_fractioninitvsfinal[i][j] = 0.0;
      avgcol_nrABperCell[i][j] = 0.0;
    }
  }

  for (int i=0; i<NRMOLTYPES; i++) {
    avgcol_totnrs[i] = 0.0;
  }
  avgcol_nrzs = 0.0;
  avgcol_nrcells = 0.0;
  avgcol_nrcells0C = 0.0;
  avgcol_nrcells1C = 0.0;
  avgcol_nrcellsmC = 0.0;
  avgcol_nrcellspC = 0.0;

  for (int i=0; i<(NRHISTBINS+1); i++) {
    avgcol_cellsizehist[i] = 0.0;
    avgcol_cellsizehist0A[i] = 0.0;
    avgcol_cellsizehist1A[i] = 0.0;
    avgcol_cellsizehistmA[i] = 0.0;
    avgcol_cellsizehistpA[i] = 0.0;
    for (int j=0; j<NRMOLTYPES; j++) {
      avgcol_nrspertypecellhist[j][i] = 0.0;
    }
  }

  return;
}


void Simbox::UpdateAveragesCollection(int, double) {
  nrhistupdates++;

  for (int i=0; i<NRMOLTYPES; i++) {
    avgcol_totnrs[i] += (double)totnrmolspertype[i];
  }
  avgcol_nrzs += (double)nr_zs;
  avgcol_nrcells += (double)nrprotocells; 

  for (int i=0; i<nrprotocells; i++) {
    int m = protocells[usednrs[i]]->GetTotalNumberOfMolecules();
    int mh = m;
    if (mh >= NRHISTBINS) { mh = NRHISTBINS; }
    avgcol_cellsizehist[mh] += 1;
    int nrB = protocells[usednrs[i]]->GetNumberOfMoleculesOfType(1);
    int nrnonB = m - nrB;

    //added by sam 
    int adjusted_nonB=nrnonB;
    //int adjusted_m=m;
    int adjusted_B=m-nrnonB;

    double myfract = (100.0*nrnonB)/m;
    avgcol_fractionpresent[(int)myfract] += 1;
    
    if (adjusted_nonB>100){
        adjusted_nonB=100;
      } 
    if (adjusted_B>100){
        adjusted_B=100;
    }
    avgcol_nrABperCell[adjusted_nonB][adjusted_B]+=1;
  
    mh = m;
    if (mh >= NRHISTBINS) { mh = NRHISTBINS; }
    if (nrnonB == 0) { 
      avgcol_cellsizehist0A[mh] += 1; 
      avgcol_nrcells0C += 1;
    }
    else if (nrnonB == 1) { 
      avgcol_cellsizehist1A[mh] += 1; 
      avgcol_nrcells1C += 1;
    }
    else { 
      avgcol_cellsizehistmA[mh] += 1; 
      avgcol_nrcellsmC += 1;
    }
    if (nrB == 0) { 
      avgcol_cellsizehistpA[mh] += 1; 
      avgcol_nrcellspC += 1;
    }
    for (int j=0; j<NRMOLTYPES; j++) {
      m = protocells[usednrs[i]]->GetNumberOfMoleculesOfType(j);
      if (m >= NRHISTBINS) { m = NRHISTBINS; }
      avgcol_nrspertypecellhist[j][m] += 1;
    }
  }

  return;
}


void Simbox::ReportAveragesCollection(void) {
  fstream of;
  double tval = 1.0/((double)nrhistupdates);


  of.open("averages.pco", ios::out);
  if(!of.is_open()) {
    throw new Exception("Could not open file 'averages.cmo'.");
  }

  for (int i=0; i<NRMOLTYPES; i++) {
    of << "<N_" << moltypenames[i] << ">:\t" << avgcol_totnrs[i]*tval << endl;
  }
  of << "<N_Z>:\t" << avgcol_nrzs*tval << endl;
  of << "<nrcells>:\t" <<  avgcol_nrcells*tval << endl;

  of << "<nrcells0C>:\t" << avgcol_nrcells0C*tval << endl;
  of << "<nrcells1C>:\t" << avgcol_nrcells1C*tval << endl;
  of << "<nrcellsmultipleC>:\t" << avgcol_nrcellsmC*tval << endl;
  of << "<nrcellspureC>:\t" << avgcol_nrcellspC*tval << endl;

  of.close();


  of.open("cellcontents.pco", ios::out);
  if(!of.is_open()) {
    throw new Exception("Could not open file 'cellcontents.cmo'.");
  }

  for (int i=0; i<101; i++) {
    for (int j=0; j<101; j++){
      of <<  avgcol_nrABperCell[i][j]*tval <<"\t";
    } 
    of << endl;
  }

  of.close();


  of.open("cellsizehist.pco", ios::out);
  if(!of.is_open()) {
    throw new Exception("Could not open file 'cellsizehist.cmo'.");
  }
  for (int i=0; i<(NRHISTBINS+1); i++) {
    of << i << "\t" << avgcol_cellsizehist[i]*tval;
    of << "\t" << avgcol_cellsizehist0A[i]*tval;
    of << "\t" << avgcol_cellsizehist1A[i]*tval;
    of << "\t" << avgcol_cellsizehistmA[i]*tval;
    for (int j=0; j<NRMOLTYPES; j++) {
      of << "\t" << avgcol_nrspertypecellhist[j][i]*tval;
    }
    of << "\t" << avgcol_cellsizehistpA[i]*tval;
    of << endl;
  }
  of.close();

  of.open("fractions.pco", ios::out);
  if(!of.is_open()) {
    throw new Exception("Could not open file 'fractions.pco'.");
  }
  double t_p = 0.0;
  double t_b = 0.0;
  double t_a = 0.0;
  for (int i=0; i<101; i++) {
    t_p += avgcol_fractionpresent[i];
    t_b += avgcol_fractionbeforedivision[i];
    t_a += avgcol_fractionafterdivision[i];
  }
  for (int i=0; i<101; i++) {
    of << i;
    of << "\t" << avgcol_fractionpresent[i]/t_p;
    of << "\t" << avgcol_fractionbeforedivision[i]/t_b;
    of << "\t" << avgcol_fractionafterdivision[i]/t_a;
    of << "\n";
  }
  of.close();


  of.open("fractioninitvsfinal.pco", ios::out);
  if(!of.is_open()) {
    throw new Exception("Could not open file 'fractioninitvsfinal.pco'.");
  }
  double t = 0.0;
  for (int i=0; i<101; i++) {
    for (int j=0; j<101; j++) {
      t += avgcol_fractioninitvsfinal[i][j];
    }
  }
  for (int i=0; i<101; i++) {
    for (int j=0; j<101; j++) {
      of << i << "\t" << j << "\t" << avgcol_fractioninitvsfinal[i][j]/t << endl;
    }
  }
  of.close();


  return;
}
