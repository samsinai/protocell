#ifndef SIMBOX
#define SIMBOX

#include <vector>
#include "protocell.h"
#include "mytypes.h"

typedef vector<Protocell*> ProtocellPtrVector;
typedef vector<Protocell*>::iterator ProtocellPtrVectorIterator;

#define gasconstantR 0.008314472 // kJ / (mol K)

using namespace std;

class Simbox {
 public:
  Simbox();
  ~Simbox();

  /*
   * Runs the simulation...
   */
  void RunSimulation(void);  
  void ReportStatus(ostream&, long&);
  void CalculateMesoscopicRateConstants(void);
  double calculate_propensities_protocell(int);
  double get_storedtotalpropensity_protocell(int nr);
  int perform_reaction(int, double, ostream&, int);
  void AddProtocell(Protocell *, int);
  void RemoveProtocell(int);
  void DivideProtocell(double, int, ostream&, int);
  void Check_propensity_consistency(void);

  void InitAveragesCollection(void);
  void UpdateAveragesCollection(int, double);
  void ReportAveragesCollection(void);

  int seed;
  long start_it;
  long n_iter;
  long erep;
  long crep;
  int timebasedrep;
  double nexttimebasedrep;
  int verbose; 
  long avgrep;
  long avgstart;
  long avgrepnext;
  long nrhistupdates;
  int logreactions;
  int splittype;
  int intcomp;

  long nextuid;

  int nrprotocells;
  int maxnrprotocells;
  double conc;
  double volume;
  double invvol;
  double time;
  int nr_zs;
  int revive;
  int modelnr;
  double divrep;

  int nrrevivals;

  double prob_aA, prob_aC, prob_aD, prob_aE, prob_aF;
  double prob_s, prob_q, prob_phi, prob_theta, prob_z, conc_z0;
  double prob_d0, prob_d1;
  int maxprotocellsize;

  ProtocellPtrVector protocells;
  vector<double> propensities;
  vector<int> availablenrs;
  vector<int> usednrs;
  vector<int> snrinfo;
 // vector< vector<int> > nrABperCell; //sam
  double summedpropensities;
  int nrpropsperprotocell;

  long nrcreations[NRMOLTYPES];
  long nrremovals[NRMOLTYPES];
  long nrprotocelldivisions;

  int totnrmolspertype[NRMOLTYPES];

  int nrBcells;

  double avgcol_totnrs[NRMOLTYPES];
  double avgcol_nrzs;
  double avgcol_nrcells;
  double avgcol_nrcells0C;
  double avgcol_nrcells1C;
  double avgcol_nrcellsmC;
  double avgcol_nrcellspC;
  double avgcol_cellsizehist[NRHISTBINS+1];
  double avgcol_cellsizehist0A[NRHISTBINS+1];
  double avgcol_cellsizehist1A[NRHISTBINS+1];
  double avgcol_cellsizehistmA[NRHISTBINS+1];
  double avgcol_cellsizehistpA[NRHISTBINS+1];
  double avgcol_nrspertypecellhist[NRMOLTYPES][NRHISTBINS+1];
  double avgcol_fractionpresent[101];
  double avgcol_fractionbeforedivision[101];
  double avgcol_fractionafterdivision[101];

  double avgcol_fractioninitvsfinal[101][101];
 
  double avgcol_nrABperCell[101][101];
  double avgcol_nrACperCell[101][101];
  double avgcol_nrADperCell[101][101];
  double avgcol_nrAEperCell[101][101];
  double avgcol_nrAFperCell[101][101];
 
  double avgcol_nrBCperCell[101][101];
  double avgcol_nrBDperCell[101][101];
  double avgcol_nrBEperCell[101][101];
  double avgcol_nrBFperCell[101][101];
 
  double avgcol_nrCDperCell[101][101];
  double avgcol_nrCEperCell[101][101];
  double avgcol_nrCFperCell[101][101];
  
  double avgcol_nrDEperCell[101][101];
  double avgcol_nrDFperCell[101][101];
  
  double avgcol_nrEFperCell[101][101];

  double avgcol_nrnonA_AperCell[101][101];  
  double avgcol_nrnonB_BperCell[101][101];
  double avgcol_nrnonC_CperCell[101][101];
  double avgcol_nrnonD_DperCell[101][101];
  double avgcol_nrnonD_EperCell[101][101];
  double avgcol_nrnonD_FperCell[101][101];

//  double all_ABCDEF[101][101][101][101][101][101][101];

};

#endif
