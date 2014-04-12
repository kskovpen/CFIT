// CFIT class definition and includes

#ifndef CFIT_H
#define CFIT_H

#include "TROOT.h"
#include "THStack.h"
#include "Riostream.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TMinuit.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TH1D.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TDecompLU.h"
#include "TDecompQRH.h"
#include "TDecompChol.h"
#include "TDecompBK.h"

#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

extern std::vector<double> sf;
extern std::vector<double> sferr;
const double PI = 3.141592653589793238462;
extern TFile *fcov;
extern int nT;
extern int nSYS;
extern int nBINS;
extern std::vector<int> vecM;
const int vecM_MAX = 100; // the max number of bins
const int TMAX = 50; // the max number of templates
const int SYSMAX = 100; // the max number of systematics
const int HMAX = 10000; // the total max number of histograms
extern TMatrixD *covMIp;
extern double cov[vecM_MAX][vecM_MAX];
extern double covFCUR;
extern double norm1D[vecM_MAX];
extern TVectorD *norm1Dp;
extern double covFactor[SYSMAX];
extern std::string nameT[TMAX];
extern std::string nameNOM;
extern std::string nameSYS[SYSMAX];
extern std::string nameSYSUP[SYSMAX];
extern std::string nameSYSLOW[SYSMAX];

extern bool VALID;
extern double CHISQ;
extern double PAR[TMAX];
extern double ERR[TMAX];

extern std::ifstream _fcom;
extern TFile *dfile;
extern TH1D *hist[HMAX];
extern TH1D *histNOM[HMAX];
extern TH1D *histLOW[HMAX];
extern TH1D *histUP[HMAX];

extern TH1D *h_comb;
extern TH1D *h_sys_low_comb[SYSMAX];
extern TH1D *h_sys_up_comb[SYSMAX];
extern TH1D *h_combNATURAL;
extern TH1D *h_sys_low_combNATURAL[SYSMAX];
extern TH1D *h_sys_up_combNATURAL[SYSMAX];
extern TH1D *h_sys_low_combSum;
extern TH1D *h_sys_up_combSum;

extern double parg[TMAX];
extern double perg[TMAX];

extern double parval[TMAX], parerr[TMAX],
  chis[10];

extern double funcDEFM(double vdata1,
		       double vdataErr1,
		       std::vector<double> vmc1,
		       double vdata2,
		       double vdataErr2,
		       std::vector<double> vmc2,
		       double covEL,
		       double norm1,
		       double norm2,
		       double *par);

extern void fcnSysM(int &npar, double *gin, double &f, double *par, int iflag);

namespace CFIT
{
   
   class cfit
     {
	
      public:
	cfit(std::string fin);
	virtual ~cfit();

	// user methods
	void run();
	void silentMode();
	double getCHISQ();
	double getPAR(int idx);
	double getERR(int idx);
	
      protected:
	
	void totSys(TH1D *h_nom,TH1D *h_low,TH1D *h_up);
	
	void combSysLinear(std::vector<TH1D*> h_nom,
			   std::vector<TH1D*> h_sys_up,
			   std::vector<TH1D*> h_sys_low,
			   TH1D *h_comb,
			   TH1D *h_sys_low_comb,
			   TH1D *h_sys_up_comb);
	
	void combSys(TH1D *h_nom,TH1D** h_sys_low,TH1D** h_sys_up,TH1D *h_sys_low_comb,TH1D *h_sys_up_comb);

	void applySys(std::vector<double> sf,
		      std::vector<double> sferr,
		      bool write);
	
	void doFit(double *par,double *err,double *chis);
	
	void doFracSys(TH1D *hnom,TH1D *hsysDown,TH1D *hsysUp,int isys);
	
	double bbCorr(TH1D *hsysDown,TH1D *hsysUp);
	
	void applySF();
	
	TGraphAsymmErrors* makeErrorBand(TH1D* tot, TH1D* plus, TH1D* minus);
	
	void removeMCSys(TH1D *hnom,TH1D *hsysUp,TH1D *hsysLow);
	
      ClassDef(CFIT::cfit,1)
     };
}

#endif
