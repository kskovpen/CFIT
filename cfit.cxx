// CFIT class implementation

#include "cfit.h"

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

// global variables and functions
bool verb = 1;
std::vector<double> sf;
std::vector<double> sferr;
double covFactor[SYSMAX];
TFile *fcov;
int nT;
int nSYS;
int asymm = 1; // 1 - max; 0 - Barlow
bool rcov;
std::string nameT[TMAX];
std::string nameNOM;
std::string nameSYS[SYSMAX];
std::string nameSYSUP[SYSMAX];
std::string nameSYSLOW[SYSMAX];
int nBINS;
std::vector<int> vecM;
TMatrixD covM(vecM_MAX,vecM_MAX);
TMatrixD covMI(vecM_MAX,vecM_MAX);
TMatrixD *covMIp;
double cov[vecM_MAX][vecM_MAX];
double norm1D[vecM_MAX];
TVectorD *norm1Dp;
TVectorD vconv(vecM_MAX);

bool VALID = 0;
double CHISQ = -1;
double PAR[TMAX];
double ERR[TMAX];

std::ifstream _fcom;
TFile *dfile;
TH1D *h_data;
TH1D *hist[HMAX];
TH1D *histNOM[HMAX];
TH1D *histLOW[HMAX];
TH1D *histUP[HMAX];

TH1D *h_comb;
TH1D *h_sys_low_comb[SYSMAX];
TH1D *h_sys_up_comb[SYSMAX];
TH1D *h_combNATURAL;
TH1D *h_sys_low_combNATURAL[SYSMAX];
TH1D *h_sys_up_combNATURAL[SYSMAX];
TH1D *h_sys_low_combSum;
TH1D *h_sys_up_combSum;

double parg[TMAX];
double perg[TMAX];

double parval[TMAX], parerr[TMAX],
  chis[10];

double funcDEFM(double vdata1,
		double vdataErr1,
		std::vector<double> vmc1,
		double vdata2,
		double vdataErr2,
		std::vector<double> vmc2,
		double covEL,
		double norm1,
		double norm2,
		double *par);

void fcnSysM(int &npar, double *gin, double &f, double *par, int iflag);

ClassImp(CFIT::cfit)

// Default constructor
// Open file and read histograms
CFIT::cfit::cfit(std::string fin)
{
   gErrorIgnoreLevel = 2000;
   
   _fcom.open(fin.c_str());
   
   if( ! _fcom.is_open() )
     {
	std::cout << "Input file does not exist: " << fin << std::endl;
	exit(1);
     }   
   
   nSYS = 0;
   rcov = 0;
   std::string line;
   while( getline(_fcom,line) )
     {
	std::vector<std::string> spl;
	split(spl,line,boost::algorithm::is_any_of(" "));

	if( strcmp(&spl[0][0],"#") == 0 ) continue;	
	else if( strcmp(spl[0].c_str(),"[FILE]") == 0 )
	  {
	     dfile = TFile::Open(spl[1].c_str());
	  }
	else if( strcmp(spl[0].c_str(),"[DATA]") == 0 )
	  {
	     h_data = (TH1D*)dfile->Get(spl[1].c_str());
	     nBINS = h_data->GetXaxis()->GetNbins();
	  }
	else if( strcmp(spl[0].c_str(),"[COV]") == 0 )
	  {
	     if( strcmp(spl[1].c_str(),"READ") == 0 )
	       {		  
		  fcov = TFile::Open(spl[1].c_str());
		  rcov = 1;
	       }
	     else
	       {
		  fcov = NULL;
	       }	     
	  }
	else if( strcmp(spl[0].c_str(),"[TEMPLATE]") == 0 )
	  {
	     nT = spl.size()-1;
	     for(int i=1;i<nT+1;i++)
	       {
		  nameT[i-1] = spl[i];
	       }	     
	  }
	else if( strcmp(spl[0].c_str(),"[NOMINAL]") == 0 )
	  {
	     nameNOM = spl[1];
	  }	
	else if( strcmp(spl[0].c_str(),"[SYS]") == 0 )
	  {
	     nameSYS[nSYS] = spl[1];
	     nameSYSLOW[nSYS] = spl[2];
	     nameSYSUP[nSYS] = spl[3];
	     nSYS++;
	  }	
     }
   
   if( nT > TMAX )
     {
	std::cout << "Max number of templates is " << TMAX << std::endl;
	exit(1);
     }   
   if( nSYS > SYSMAX )
     {
	std::cout << "Max number of systematic variations is " << SYSMAX << std::endl;
	exit(1);
     }   
   if( nBINS > vecM_MAX )
     {
	std::cout << "Max number of bins in template is " << vecM_MAX << std::endl;
	exit(1);
     }   

   for(int it=0;it<nT;it++)
     {	
	std::string hnameNOM = nameT[it]+nameNOM;
	TH1D *htNOM = (TH1D*)dfile->Get(hnameNOM.c_str());
	std::string hcopyNOM = hnameNOM+"_NOMcopy";
	histNOM[it] = (TH1D*)htNOM->Clone(hcopyNOM.c_str());
	
	for(int is=0;is<nSYS;is++)
	  {
	     std::string hnameLOW = nameT[it]+nameSYSLOW[is];
	     TH1D *htLOW = (TH1D*)dfile->Get(hnameLOW.c_str());
	     std::string hcopyLOW = hnameLOW+"_copy";
	     histLOW[is+it*nSYS] = (TH1D*)htLOW->Clone(hcopyLOW.c_str());

	     std::string hnameUP = nameT[it]+nameSYSUP[is];
	     TH1D *htUP = (TH1D*)dfile->Get(hnameUP.c_str());
	     std::string hcopyUP = hnameUP+"_copy";
	     histUP[is+it*nSYS] = (TH1D*)htUP->Clone(hcopyUP.c_str());
	  }
	
	if( nSYS == 0 )
	  {
	     std::string hnameLOW = nameT[it]+nameNOM;
	     TH1D *htLOW = (TH1D*)dfile->Get(hnameLOW.c_str());
	     std::string hcopyLOW = hnameLOW+"_copy";
	     histLOW[it] = (TH1D*)htLOW->Clone(hcopyLOW.c_str());

	     std::string hnameUP = nameT[it]+nameNOM;
	     TH1D *htUP = (TH1D*)dfile->Get(hnameUP.c_str());
	     std::string hcopyUP = hnameUP+"_copy";
	     histUP[it] = (TH1D*)htUP->Clone(hcopyUP.c_str());
	  }	
     }
}

// Destructor
CFIT::cfit::~cfit()
{
   fcov->Close();
}

// chi2 term calculation
double funcDEFM(double vdata1,
			    double vdataErr1,
			    std::vector<double> vmc1,
			    double vdata2,
			    double vdataErr2,
			    std::vector<double> vmc2,
			    double covEL,
			    double norm1,
			    double norm2,
			    double *par)
{
   double term1 = 0.;
   double term2 = 0.;
   
   for(int i=0;i<vmc1.size();i++)
     {
	term1 += par[i]*vmc1[i];
	term2 += par[i]*vmc2[i];
     }
   
   double val = ((vdata1-term1)/norm1)*((vdata2-term2)/norm2)*covEL;
   
   return val;
}

// chi2 definition
void fcnSysM(int &npar, double *gin, double &f, double *par, int iflag)
{
   double chisq = 0;
   double delta;

   double nb = vecM.size();
   
   for( int i1=0;i1<nb;i1++ )
     {
	int idx1 = vecM[i1];
	
	std::vector<double> vmc1;
	vmc1.clear();
	for(int it1=0;it1<nT;it1++)
	  {
	     vmc1.push_back(histNOM[it1]->GetBinContent(idx1));
	  }	
	
	double data_1 = h_data->GetBinContent(idx1);
	double data_err_1 = h_data->GetBinError(idx1);
	
	for( int i2=0;i2<nb;i2++ )
	  {
	     int idx2 = vecM[i2];

	     std::vector<double> vmc2;
	     vmc2.clear();
	     for(int it2=0;it2<nT;it2++)
	       {
		  vmc2.push_back(histNOM[it2]->GetBinContent(idx2));
	       }		     
	     
	     double data_2 = h_data->GetBinContent(idx2);
	     double data_err_2 = h_data->GetBinError(idx2);

	     delta = funcDEFM(data_1,
			      data_err_1,
			      vmc1,
			      data_2,
			      data_err_2,
			      vmc2,
			      (*covMIp)[i1][i2],
			      (*norm1Dp)[i1],
			      (*norm1Dp)[i2],
			      par);

	     chisq += delta;
	  }
     }
   
//   if( iflag == 1 )
//     {		
//	if( verb ) std::cout << "CHISQ/NDOF (NO FIT) = "<< chisq/(nBINS-nT) << std::endl;
//     }
   if( iflag == 3 )
     {
	VALID = 1;
	CHISQ = chisq/(nBINS-nT);
	if( verb ) std::cout << "CHISQ/NDOF = " << CHISQ << std::endl;
     }
   
   f = chisq;
}

// measure max and min systematic variations
void CFIT::cfit::totSys(TH1D *h_nom,TH1D *h_low,TH1D *h_up)
{   
   int nbins = h_nom->GetXaxis()->GetNbins();
   
   for(int ib=1;ib<=nbins;ib++)
     {
	double b_nom = h_nom->GetBinContent(ib);
	double b_cur_low = h_low->GetBinContent(ib);
	double b_cur_up = h_up->GetBinContent(ib);
	
	double del_low = b_cur_low - b_nom;
	double del_up = b_cur_up - b_nom;
	double del_low_res = del_low;
	double del_up_res = del_up;
	
	double del_min = (del_low < del_up) ? del_low : del_up;
	double del_max = (del_low > del_up) ? del_low : del_up;

	del_low_res = (del_min > 0.) ? 0. : del_min;
	del_up_res = (del_max < 0.) ? 0. : del_max;
	
	double sys_low = b_nom + del_low_res;
	double sys_up = b_nom + del_up_res;
	
	h_low->SetBinContent(ib,sys_low);
	h_up->SetBinContent(ib,sys_up);
     }
}

// quadratic combination of systematics 
void CFIT::cfit::combSys(TH1D *h_nom,TH1D** h_sys_low,TH1D** h_sys_up,TH1D *h_sys_low_comb,TH1D *h_sys_up_comb)
{
   double sys_up[10000];
   double sys_low[10000];
   
   for(int is=1;is<=nBINS;is++)
     {
	sys_up[is-1] = 0.;
	sys_low[is-1] = 0.;
     }      
   
   for(int j=0;j<nSYS;j++)
     {
	for(int ib=1;ib<=nBINS;ib++)
	  {
	     double nom = h_nom->GetBinContent(ib);
	     double low = h_sys_low[j]->GetBinContent(ib);
	     double up = h_sys_up[j]->GetBinContent(ib);
	     double delta_low = nom - low;
	     double delta_up = up - nom;
	     sys_up[ib-1] = sqrt(pow(sys_up[ib-1],2) + pow(delta_up,2));
	     sys_low[ib-1] = sqrt(pow(sys_low[ib-1],2) + pow(delta_low,2));
	  }
     }      
   
   for(int ib=1;ib<=nBINS;ib++)
     {
	double nom = h_nom->GetBinContent(ib);
	double err = h_nom->GetBinError(ib);
	
	sys_up[ib-1] = sqrt(pow(sys_up[ib-1],2) + pow(err,2));
	sys_low[ib-1] = sqrt(pow(sys_low[ib-1],2) + pow(err,2));
     }   
   
   for(int ib=1;ib<=nBINS;ib++)
     {		
	double nom = h_nom->GetBinContent(ib);
	h_sys_low_comb->SetBinContent(ib,nom-sys_low[ib-1]);
	h_sys_up_comb->SetBinContent(ib,nom+sys_up[ib-1]);
     }
}

// linear combination of systematics
void CFIT::cfit::combSysLinear(std::vector<TH1D*> h_nom,
			       std::vector<TH1D*> h_sys_up,
			       std::vector<TH1D*> h_sys_low,
			       TH1D *h_comb,
			       TH1D *h_sys_low_comb,
			       TH1D *h_sys_up_comb)
{
   double sys_up[10000];
   double sys_low[10000];
   double stat[10000];
   
   for(int is=1;is<=nBINS;is++)
     {
	sys_up[is-1] = 0.;
	sys_low[is-1] = 0.;
	stat[is-1] = 0.;
     }

   for(int ib=1;ib<=nBINS;ib++)
     {
	double delta_low = 0.;
	double delta_up = 0.;
	
	for(int i=0;i<h_nom.size();i++)
	  {
	     double nom = h_nom[i]->GetBinContent(ib);
	     double up = h_sys_up[i]->GetBinContent(ib);
	     double low = h_sys_low[i]->GetBinContent(ib);
	     double err = h_nom[i]->GetBinError(ib);
	     
	     delta_low += (nom-low);
	     delta_up += (up-nom);
	     
	     stat[ib-1] += pow(err,2);
	  }	
	
	sys_up[ib-1] = sys_up[ib-1] + delta_up;
	sys_low[ib-1] = sys_low[ib-1] + delta_low;

	stat[ib-1] = sqrt(stat[ib-1]);
     }      
   
   for(int ib=1;ib<=nBINS;ib++)
     {
	double nom = 0.;
	for(int i=0;i<h_nom.size();i++)
	  {
	     nom += h_nom[i]->GetBinContent(ib);
	  }	

	h_comb->SetBinContent(ib,nom);
	h_comb->SetBinError(ib,stat[ib-1]);
	h_sys_low_comb->SetBinContent(ib,nom-sys_low[ib-1]);
	h_sys_up_comb->SetBinContent(ib,nom+sys_up[ib-1]);
     }
}

// analysis of systematic uncertainties and correlation matrix construction
void CFIT::cfit::applySys(std::vector<double> sf,
			  std::vector<double> sferr,
			  bool write)
{
   std::vector<TH1D*> v_fit;
   
   for(int i=0;i<nT;i++)
     {
	std::string c_fit = nameT[i] + "_fit";
	TH1D *h_fit = (TH1D*)histNOM[i]->Clone(c_fit.c_str());
	h_fit->Scale(sf[i]);
	v_fit.push_back(h_fit);
     }   

   std::vector<std::vector<TH1D*> > v_sys_low_fit;
   std::vector<std::vector<TH1D*> > v_sys_up_fit;

   for(int i=0;i<nSYS;i++)
     {
	std::string c_sys_up_fit = nameSYS[i] + "_up_fit";
	std::string c_sys_low_fit = nameSYS[i] + "_low_fit";
   
	std::vector<TH1D*> v_sys_low;
	std::vector<TH1D*> v_sys_up;
	
	for(int j=0;j<nT;j++)
	  {
	     TH1D *h_low = (TH1D*)histLOW[i+nSYS*j]->Clone((nameT[j]+c_sys_low_fit).c_str());
	     TH1D *h_up = (TH1D*)histUP[i+nSYS*j]->Clone((nameT[j]+c_sys_up_fit).c_str());
	
	     h_low->Scale(sf[j]);
	     h_up->Scale(sf[j]);
	     
	     v_sys_low.push_back(h_low);
	     v_sys_up.push_back(h_up);
	  }
	
	v_sys_low_fit.push_back(v_sys_low);
	v_sys_up_fit.push_back(v_sys_up);
     }

   if( nSYS == 0 )
     {
	std::string c_sys_up_fit = "nosys_up_fit";
	std::string c_sys_low_fit = "nosys_low_fit";
   
	std::vector<TH1D*> v_sys_low;
	std::vector<TH1D*> v_sys_up;
	
	for(int j=0;j<nT;j++)
	  {
	     TH1D *h_low = (TH1D*)histLOW[j]->Clone((nameT[j]+c_sys_low_fit).c_str());
	     TH1D *h_up = (TH1D*)histUP[j]->Clone((nameT[j]+c_sys_up_fit).c_str());
	
	     h_low->Scale(sf[j]);
	     h_up->Scale(sf[j]);
	     
	     v_sys_low.push_back(h_low);
	     v_sys_up.push_back(h_up);
	  }
	
	v_sys_low_fit.push_back(v_sys_low);
	v_sys_up_fit.push_back(v_sys_up);
     }
      
   if( !h_comb ) delete h_comb;
   h_comb = (TH1D*)histNOM[0]->Clone("h_comb");
   if( !h_sys_low_combSum ) delete h_sys_low_combSum;
   h_sys_low_combSum = (TH1D*)histNOM[0]->Clone("h_sys_low_combSum");
   if( !h_sys_up_combSum ) delete h_sys_up_combSum;
   h_sys_up_combSum = (TH1D*)histNOM[0]->Clone("h_sys_up_combSum");
   
   for(int j=0;j<nSYS;j++)
     {
	std::string hnameLOW = "h_sys_low_comb_"+nameSYS[j];
	if( !h_sys_low_comb[j] ) delete h_sys_low_comb[j];
	h_sys_low_comb[j] = (TH1D*)h_comb->Clone(hnameLOW.c_str());

	std::string hnameUP = "h_sys_up_comb_"+nameSYS[j];
	if( !h_sys_up_comb[j] ) delete h_sys_up_comb[j];
	h_sys_up_comb[j] = (TH1D*)h_comb->Clone(hnameUP.c_str());
     }   

   if( nSYS == 0 )
     {
	std::string hnameLOW = "h_sys_low_comb";
	if( !h_sys_low_comb[0] ) delete h_sys_low_comb[0];
	h_sys_low_comb[0] = (TH1D*)h_comb->Clone(hnameLOW.c_str());

	std::string hnameUP = "h_sys_up_comb";
	if( !h_sys_up_comb[0] ) delete h_sys_up_comb[0];
	h_sys_up_comb[0] = (TH1D*)h_comb->Clone(hnameUP.c_str());
     }   
   
   for(int j=0;j<nSYS;j++)
     {
	combSysLinear(v_fit,
		      v_sys_up_fit[j],
		      v_sys_low_fit[j],
		      h_comb,
		      h_sys_low_comb[j],
		      h_sys_up_comb[j]);
	
	std::string nat_low = "h_sys_low_combNATURAL_"+nameSYSLOW[j];
	if( !h_sys_low_combNATURAL[j] ) delete h_sys_low_combNATURAL[j];
	h_sys_low_combNATURAL[j] = (TH1D*)h_sys_low_comb[j]->Clone(nat_low.c_str());
	std::string nat_up = "h_sys_up_combNATURAL_"+nameSYSUP[j];
	if( !h_sys_up_combNATURAL[j] ) delete h_sys_up_combNATURAL[j];
	h_sys_up_combNATURAL[j] = (TH1D*)h_sys_up_comb[j]->Clone(nat_up.c_str());

	if( !h_combNATURAL ) delete h_combNATURAL;
	h_combNATURAL = (TH1D*)h_comb->Clone("h_combNATURAL");

	doFracSys(h_comb,h_sys_low_combNATURAL[j],h_sys_up_combNATURAL[j],j);
	
	totSys(h_comb,
	       h_sys_low_comb[j],
	       h_sys_up_comb[j]
	      );
	
	removeMCSys(h_comb,
		    h_sys_up_comb[j],
		    h_sys_low_comb[j]);
     }

   if( nSYS == 0 )
     {
	combSysLinear(v_fit,
		      v_sys_up_fit[0],
		      v_sys_low_fit[0],
		      h_comb,
		      h_sys_low_comb[0],
		      h_sys_up_comb[0]);

	std::string nat_low = "h_sys_low_combNATURAL";
	if( !h_sys_low_combNATURAL[0] ) delete h_sys_low_combNATURAL[0];
	h_sys_low_combNATURAL[0] = (TH1D*)h_sys_low_comb[0]->Clone(nat_low.c_str());
	std::string nat_up = "h_sys_up_combNATURAL";
	if( !h_sys_up_combNATURAL[0] ) delete h_sys_up_combNATURAL[0];
	h_sys_up_combNATURAL[0] = (TH1D*)h_sys_up_comb[0]->Clone(nat_up.c_str());

	if( !h_combNATURAL ) delete h_combNATURAL;
	h_combNATURAL = (TH1D*)h_comb->Clone("h_combNATURAL");

	doFracSys(h_comb,h_sys_low_combNATURAL[0],h_sys_up_combNATURAL[0],0);
	
	totSys(h_comb,
	       h_sys_low_comb[0],
	       h_sys_up_comb[0]
	      );
	
	removeMCSys(h_comb,
		    h_sys_up_comb[0],
		    h_sys_low_comb[0]);
     }
   
   combSys(h_comb,
	   h_sys_low_comb,
	   h_sys_up_comb,
	   h_sys_low_combSum,
	   h_sys_up_combSum);
   
   totSys(h_comb,
	  h_sys_low_combSum,
	  h_sys_up_combSum);

   removeMCSys(h_comb,
	       h_sys_up_combSum,
	       h_sys_low_combSum);
   
   vecM.clear();
   
   for(int i=1;i<=nBINS;i++)
     {
	double cont = h_comb->GetBinContent(i);
	double contUp = h_sys_up_combSum->GetBinContent(i);
	double contDown = h_sys_low_combSum->GetBinContent(i);
	
	vecM.push_back(i);
     } 
   
   if( ! rcov )
     {	
	const int vecM_n = vecM.size();
	
	double norm[vecM_MAX][vecM_MAX];
	
	for(int i1=0;i1<vecM_MAX;i1++)
	  {
	     for(int i2=0;i2<vecM_MAX;i2++)
	       {
		  norm[i1][i2] = 0.;
	       }	
	  }   

	for(int i1=0;i1<vecM_MAX;i1++)
	  {
	     for(int i2=0;i2<vecM_MAX;i2++)
	       {
		  cov[i1][i2] = 0.;
	       }
	  }

	// compute covariance matrix
	for(int i1=0;i1<vecM_n;i1++)
	  {
	     int idx1 = vecM[i1];
	     
	     double siglin = 0.;
	     
	     for(int i2=0;i2<vecM_n;i2++)
	       {
		  int idx2 = vecM[i2];
		  
		  double cov_v = 0.;
		  double norm_i1_v = 0.;
		  double norm_i2_v = 0.;
		  
		  for(int is=-1;is<nSYS;is++)
		    {
		       TH1D *nom = h_combNATURAL;
		       TH1D *low = h_combNATURAL;
		       TH1D *up = h_combNATURAL;
		       TH1D *lowSum = h_sys_low_combSum;
		       TH1D *upSum = h_sys_up_combSum;

		       if( is != -1 )
			 {
			    low = h_sys_low_combNATURAL[is];
			    up = h_sys_up_combNATURAL[is];
			 }
		       
		       if( nSYS == 0 )
			 {
			    low = h_sys_low_combNATURAL[0];
			    up = h_sys_up_combNATURAL[0];
			 }		       
		       
		       double v_nom_i1 = nom->GetBinContent(idx1);
		       double v_sysUp_i1 = up->GetBinContent(idx1);
		       double v_sysDown_i1 = low->GetBinContent(idx1);
		       
		       double v_nom_i2 = nom->GetBinContent(idx2);
		       double v_sysUp_i2 = up->GetBinContent(idx2);
		       double v_sysDown_i2 = low->GetBinContent(idx2);
		       
		       double sysUp_i1 = v_sysUp_i1 - v_nom_i1;
		       double sysDown_i1 = v_nom_i1 - v_sysDown_i1;
		       
		       double sysUp_i2 = v_sysUp_i2 - v_nom_i2;
		       double sysDown_i2 = v_nom_i2 - v_sysDown_i2;
		       
		       if( asymm == 0 )
			 {			    
			    if( sysUp_i1*sysDown_i1 < 0 )
			      {
				 if( fabs(sysUp_i1) > fabs(sysDown_i1) )
				   sysDown_i1 = 0.;
				 else
				   sysUp_i1 = 0.;
			      }		  
			    if( sysUp_i2*sysDown_i2 < 0 )
			      {
				 if( fabs(sysUp_i2) > fabs(sysDown_i2) )
				   sysDown_i2 = 0.;
				 else
				   sysUp_i2 = 0.;
			      }		  
			 }		       
		       
		       if( i2 == 0 )
			 siglin += sysUp_i1 - sysDown_i1;
		       
		       double stat_i1 = nom->GetBinError(idx1);
		       
		       double stat_i2 = nom->GetBinError(idx2);
		       
		       double data_stat = h_data->GetBinError(idx1);
		       if( i1 != i2 || is != -1 ) data_stat = 0.;
		       
		       double sigma_i1 = (fabs(sysUp_i1) > fabs(sysDown_i1)) ? sysUp_i1 : sysDown_i1;
		       double sigma_i2 = (fabs(sysUp_i2) > fabs(sysDown_i2)) ? sysUp_i2 : sysDown_i2;

		       if( asymm == 0 )
			 {
			    sigma_i1 = sysUp_i1+sysDown_i1;
			    sigma_i2 = sysUp_i2+sysDown_i2;
			 }		       
		       
		       double alpha_i1 = sysUp_i1-sysDown_i1;
		       double alpha_i2 = sysUp_i2-sysDown_i2;
		       
		       sysUp_i1 = upSum->GetBinContent(idx1) - 
			 h_comb->GetBinContent(idx1);
		       sysDown_i1 = -lowSum->GetBinContent(idx1) +
			 h_comb->GetBinContent(idx1);
		       
		       sysUp_i2 = upSum->GetBinContent(idx2) - 
			 h_comb->GetBinContent(idx2);
		       sysDown_i2 = -lowSum->GetBinContent(idx2) +
			 h_comb->GetBinContent(idx2);
		       
		       double sigmaAv_i1 = (sysUp_i1+sysDown_i1)/2.;
		       double sigmaAv_i2 = (sysUp_i2+sysDown_i2)/2.;
		       
		       double covAS = 0.;
		       if( is != -1 ) covAS = covFactor[is];
		       
		       if( i1 == i2 ) covAS = 1.;
		       
		       if( asymm == 0 )
			 {		       
			    cov_v += covAS*sigma_i1/2.*sigma_i2/2.+
			      covAS*alpha_i1/2.*alpha_i2/2.*(1.-2./PI)+
			      pow(data_stat,2);
			 }
		       else if( asymm == 1 )
			 {		       
			    cov_v += covAS*sigma_i1*sigma_i2+
			      pow(data_stat,2);
			 }		  
		       
		       if( i1 == i2 && is == -1 )
			 cov_v += pow(stat_i1,2);
		       
		       if( is != -1 )
			 {		       
			    stat_i1 = 0;
			    stat_i2 = 0;
			 }
		       
		       if( asymm == 0 )
			 {		       
			    norm_i1_v += covAS*sigma_i1/2.*sigma_i1/2.+
			      covAS*alpha_i1/2.*alpha_i1/2.*(1.-2./PI) +
			      pow(stat_i1,2)+pow(data_stat,2);
			    
			    norm_i2_v += covAS*sigma_i2/2.*sigma_i2/2.+
			      covAS*alpha_i2/2.*alpha_i2/2.*(1.-2./PI) +
			      pow(stat_i2,2)+pow(data_stat,2);
			 }
		       else if( asymm == 1 )
			 {		       
			    norm_i1_v += covAS*sigma_i1*sigma_i1+
			      pow(stat_i1,2)+pow(data_stat,2);
			    norm_i2_v += covAS*sigma_i2*sigma_i2+
			      pow(stat_i2,2)+pow(data_stat,2);
			 }
		    }	     
		  
		  norm[i1][i2] = sqrt(norm_i1_v)*sqrt(norm_i2_v);
		  
		  if( i1 == i2 ) 
		    {
		       norm1D[i1] = sqrt(norm_i1_v);
		    }

		  cov[i1][i2] = cov_v;
	       }
	  }
	
	for(int i1=0;i1<vecM_n;i1++)
	  {
	     for(int i2=0;i2<vecM_n;i2++)
	       {
		  cov[i1][i2] = cov[i1][i2] / norm1D[i1] / norm1D[i2];		  
	       }
	  }   		
	
	covM.ResizeTo(vecM_n,vecM_n);
	covMI.ResizeTo(vecM_n,vecM_n);
	for(int i1=0;i1<vecM_n;i1++)
	  {
	     for(int i2=0;i2<vecM_n;i2++)
	       {
		  covM[i1][i2] = cov[i1][i2];
	       }	
	  }
	
	TMatrixDSym covMSym;
	TMatrixDSym covMISym;
	covMSym.ResizeTo(vecM_n,vecM_n);
	covMISym.ResizeTo(vecM_n,vecM_n);
	for(int i1=0;i1<vecM_n;i1++)
	  {
	     for(int i2=0;i2<vecM_n;i2++)
	       {
		  covMSym[i1][i2] = cov[i1][i2];
		  covMISym[i1][i2] = cov[i1][i2];
	       }	
	  }	
	TDecompBK lu(covMSym);
	lu.Invert(covMISym);

	TMatrixD unityM = covMSym*covMISym;

	if( fabs(unityM.Determinant()-1) > 10E-6 )
	  {
	     std::cout << "Inversion of covariance matrix failed" << std::endl;
	     std::cout << unityM.Determinant()-1 << std::endl;
	     exit(1);
	  }      

	for(int i1=0;i1<vecM_n;i1++)
	  {
	     for(int i2=0;i2<vecM_n;i2++)
	       {
		  covMI[i1][i2] = covMISym[i1][i2];
	       }	
	  }
	
	if( covMI.Determinant() < 0 )
	  {
	     std::cout << "Determinant is negative" << std::endl;
	     exit(1);
	  }

	for(int ie=0;ie<vecM_MAX;ie++)
	  {
	     vconv[ie] = norm1D[ie];
	  }

	norm1Dp = &vconv;

	if( write )
	  {	     
	     fcov = new TFile("cov.root","RECREATE");

	     std::string matrixName = "covMI";
	     covMI.Write(matrixName.c_str());
	     std::string normName = "norm";
	     vconv.Write(normName.c_str());
	  }
	     
	covMIp = &covMI;
     }
   else
     {
	std::string matrixName = "covMI";
	covMIp = (TMatrixD*)fcov->Get(matrixName.c_str());
	std::string normName = "norm";
	norm1Dp = (TVectorD*)fcov->Get(normName.c_str());
     }

   for(int i=0;i<nT;i++)
     {
	delete v_fit[i];
     }
   
   for(int i=0;i<nSYS;i++)
     {
	for(int j=0;j<nT;j++)
	  {
	     delete v_sys_low_fit[i][j];
	     delete v_sys_up_fit[i][j];
	  }	
     }
}

// minuit fit
void CFIT::cfit::doFit(double *par,double *err,double *chis)
{
   double pval[nT], perr[nT];
   
   TMinuit *gMinuit = new TMinuit(1);
   gMinuit->SetFCN(fcnSysM);
   gMinuit->SetPrintLevel(-1);
   
   Double_t arglist[10];
   Int_t ierflg = 0;
   
   arglist[0] = 3;
   gMinuit->mnexcm("SET ERR",arglist,1,ierflg);   
   
   for(int i=0;i<nT;i++)
     {
	std::string pname = "p"+std::string(Form("%d",i));
	gMinuit->mnparm(i,pname.c_str(),1.,pow(10.,-4),0.00001,100,ierflg);
	parg[i] = 1.;
	perg[i] = 0.;
     }

/*   
   if( h_draw_cjet[iVARcur]->Integral() < 0.001 )
     {
	if( h_draw_cjet[iVARcur]->GetEntries() < 3 )
	  gMinuit->mnparm(1,"p2",0.,step[1],0.00001,100,ierflg);

	gMinuit->FixParameter(1);
     }   
   if( h_draw_ljet[iVARcur]->Integral() < 0.001 )
     {
	if( h_draw_ljet[iVARcur]->GetEntries() < 3 )
	  gMinuit->mnparm(2,"p3",0.,step[2],0.00001,100,ierflg);
	
	gMinuit->FixParameter(2);
     }   
*/  
   gMinuit->mnexcm("MIGRAD", arglist ,0,ierflg);

   for(int i=0;i<nT;i++)
     {   
	gMinuit->GetParameter(i,parg[i],perg[i]);
//	std::cout << parg[i] << "+" << perg[i] << std::endl;
     }   

   sf.clear();
   sferr.clear();
   for(int i=0;i<nT;i++)
     {
	sf.push_back(parg[i]);
	sferr.push_back(0);
     }   
   
   applySys(sf,sferr,1);
   
   gMinuit->mnexcm("MIGRAD", arglist ,0,ierflg);
   gMinuit->mnexcm("MINOS", arglist ,0,ierflg);
   gMinuit->mnexcm("MIGRAD", arglist ,0,ierflg);
   gMinuit->mnexcm("MINOS", arglist ,0,ierflg);

   double paramErrP[nT];
   double paramErrM[nT];
   double paramErrS[nT];
   double paramErrG[nT];
  
   for(int i=0;i<nT;i++)
     {	
	gMinuit->mnerrs(i,paramErrP[i],paramErrM[i],paramErrS[i],paramErrG[i]);
     }   
   
   arglist[0] = 3;
   gMinuit->mnexcm("CALL FCN",arglist,1,ierflg);
   
   double amin,edm,errdef;
   int nvpar,nparx,icstat;   
   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
   chis[0] = amin;
   chis[1] = nBINS - gMinuit->GetNumFreePars();

   for(int i=0;i<nT;i++)
     {	   
	gMinuit->GetParameter(i,pval[i],perr[i]);
	par[i] = pval[i];
	err[i] = perr[i];
	if( verb ) std::cout << nameT[i] << " " << par[i] << " +- " << err[i] << std::endl;
     }   

   sf.clear();
   sferr.clear();
   for(int i=0;i<nT;i++)
     {
	sf.push_back(par[i]);
	sferr.push_back(err[i]);
     }   
   
   delete gMinuit;
}

// user main function to perform a fit
void CFIT::cfit::run()
{
   sf.clear();
   sferr.clear();
   
   for(int i=0;i<nT;i++)
     {
	sf.push_back(1.);
	sferr.push_back(0.);
     }
   
   applySys(sf,sferr,0);

   doFit(parval,parerr,chis);

   applySys(sf,sferr,0);
   
   applySF();
}

// analysis of fractional systematic uncertainties
void CFIT::cfit::doFracSys(TH1D *hnom,TH1D *hsysDown,TH1D *hsysUp,int isys)
{
   TH1D *hsysDownf = (TH1D*)hnom->Clone("hsysDownf");
   TH1D *hsysUpf = (TH1D*)hnom->Clone("hsysUpf");

   TH1D *hsysDownf_cov = (TH1D*)hnom->Clone("hsysDownf_cov");
   TH1D *hsysUpf_cov = (TH1D*)hnom->Clone("hsysUpf_cov");
   
   for(int ib0=1;ib0<nBINS+1;ib0++)
     {
	float b0 = hnom->GetBinContent(ib0);
	float b0err = hnom->GetBinError(ib0);
	float bDown = hsysDown->GetBinContent(ib0);
	float bUp = hsysUp->GetBinContent(ib0);
	float bDown_cov = bDown;
	float bUp_cov = bUp;
	
	if( b0 > 0. )
	  {	     
	     bDown = (b0-bDown)/b0;
	     bUp = (b0-bUp)/b0;
	     
	     bDown_cov = (b0-bDown_cov)/b0;
	     bUp_cov = (bUp_cov-b0)/b0;
	  }
	else
	  {
	     bDown = 0.;
	     bUp = 0.;

	     bDown_cov = 0.;
	     bUp_cov = 0.;
	  }
	
	hsysDownf->SetBinContent(ib0,1.-bDown);
	hsysDownf->SetBinError(ib0,0.);
	hsysUpf->SetBinContent(ib0,1.-bUp);
	hsysUpf->SetBinError(ib0,0.);

	hsysDownf_cov->SetBinContent(ib0,bDown_cov);
	hsysDownf_cov->SetBinError(ib0,0.);
	hsysUpf_cov->SetBinContent(ib0,bUp_cov);
	hsysUpf_cov->SetBinError(ib0,0.);
     }	

   double corr = bbCorr(hsysDownf_cov,hsysUpf_cov);
   covFactor[isys] = corr;

   TCanvas *c1 = new TCanvas();
   gStyle->SetOptStat(0);
   c1->Clear();
   c1->SetLogy(0);
   hsysDownf->SetFillStyle(0);
   hsysUpf->SetFillStyle(0);
   hsysDownf->SetLineColor(kBlue);
   hsysDownf->SetLineStyle(9);
   hsysDownf->SetLineWidth(2);
   hsysUpf->SetLineColor(kRed);
   hsysUpf->SetLineStyle(1);
   hsysUpf->SetLineWidth(2);
   hsysDownf->Draw("hist");
   hsysUpf->Draw("hist same");
   double xmin = hsysUpf->GetXaxis()->GetBinUpEdge(nBINS);
   double xmax = hsysUpf->GetXaxis()->GetBinLowEdge(1);
   TLine *l1 = new TLine(xmin,1.0,xmax,1.0);
   l1->SetLineStyle(2);
   l1->SetLineWidth(1);
   l1->Draw("same");

   hsysDownf->GetXaxis()->SetTitle("Fitting variable");
   hsysDownf->GetYaxis()->SetTitle("#sigma(A)/A");
   hsysDownf->SetTitle(nameSYS[isys].c_str());      
   
   TLegend *legf = new TLegend(0.70,0.85,0.85,0.68);
   legf->SetFillColor(253);
   legf->SetBorderSize(0);   
   legf->AddEntry(hsysDownf,"- #sigma","l");
   legf->AddEntry(hsysUpf,"+ #sigma","l");
   legf->Draw();
   
   double ymax = 0.;
   double ymin = 10.;
   for(int i=1;i<=nBINS;i++)
     {
	double y1 = hsysDownf->GetBinContent(i);
	double y2 = hsysUpf->GetBinContent(i);
	
	if( y1 < ymin ) ymin = y1;	
	if( y2 < ymin ) ymin = y2;

	if( y1 > ymax ) ymax = y1;
	if( y2 > ymax ) ymax = y2;
     }   

   hsysDownf->GetYaxis()->SetRangeUser(1.-2.0*(1.-ymin),2.0*(ymax-1.)+1.);
   
   std::string picname = "SYS.eps(";
   if( isys > 0 && isys != nSYS-1 ) picname = "SYS.eps";
   else if( isys == nSYS-1 ) picname = "SYS.eps)";
   c1->Print(picname.c_str());
   c1->Clear();
   delete l1;
   delete legf;
   delete c1;
}

// computation of correlation factor for each systematics
double CFIT::cfit::bbCorr(TH1D *hsysDown,TH1D *hsysUp)
{
   std::vector<int> bb;
   bb.clear();
   
   for(int ib=1;ib<=nBINS;ib++)
     {
	double v1 = hsysDown->GetBinContent(ib);
	double v2 = hsysUp->GetBinContent(ib);
	
	if( v1 != 0 || v2 != 0 ) bb.push_back(ib);
     }   
   
   int sz = bb.size();
         
   double rho = 0.;

   if( sz < 2 ) return 1.;
   
   for(int ib1=1;ib1<=sz-1;ib1++)
     {
	for(int ib2=ib1+1;ib2<=sz;ib2++)
	  {	     
	     int id1 = bb[ib1-1];
	     int id2 = bb[ib2-1];
	     
	     double si_Down = hsysDown->GetBinContent(id1);
	     double sj_Down = hsysDown->GetBinContent(id2);
	     
	     double si_Up = hsysUp->GetBinContent(id1);
	     double sj_Up = hsysUp->GetBinContent(id2);
	     
	     double rhoij = (si_Down*sj_Down + si_Up*sj_Up);
	     double norm1 = sqrt(si_Down*si_Down+si_Up*si_Up);
	     double norm2 = sqrt(sj_Down*sj_Down+sj_Up*sj_Up);
	     if( norm1 > 0 && norm2 > 0 ) rhoij /= (norm1*norm2);
	     else rhoij = 0.;	     
	     
	     rho += fabs(rhoij);
	  }
     }   
   
   rho /= (sz*(sz-1)/2);
   
   return rho;
}

// apply fit parameters to the templates
void CFIT::cfit::applySF()
{
   TCanvas *c1 = new TCanvas();
   c1->Clear();
   
   THStack *h_draw_st_fit = new THStack();

   TLegend *legf = new TLegend(0.70,0.85,0.85,0.68);
   legf->SetFillColor(253);
   legf->SetBorderSize(0);   
   
   TH1D *h_clone[nT];
   for(int i=0;i<nT;i++)
     {	
	std::string hname = "hist"+nameT[i]+"_clone";
	h_clone[i] = (TH1D*)(histNOM[i])->Clone(hname.c_str());
	h_clone[i]->Scale(sf[i]);
	PAR[i] = sf[i];
	ERR[i] = sferr[i];
	h_clone[i]->SetMarkerSize(0);
	h_clone[i]->SetLineColor(2+i);
	h_clone[i]->SetFillColor(2+i);
	h_draw_st_fit->Add(h_clone[i]);

	legf->AddEntry(h_clone[i],nameT[i].c_str(),"f");
     }   

   h_data->SetMarkerStyle(20);
   h_data->SetMarkerSize(0.8);
   h_draw_st_fit->Draw("hist e1");
   h_data->Draw("e1 same");
   
   TGraphAsymmErrors *gr_mc_merged = makeErrorBand(h_comb,
						   h_sys_up_combSum,
						   h_sys_low_combSum);

   h_comb->SetLineColor(kBlack);
   h_comb->Draw("hist same");

   gr_mc_merged->SetFillStyle(3004);
   gStyle->SetHatchesLineWidth(2);
   gr_mc_merged->Draw("2SAME");

   h_draw_st_fit->GetXaxis()->SetTitle("Fitting variable");
   h_draw_st_fit->GetYaxis()->SetTitle("Fractions");
   
   double max = h_draw_st_fit->GetMaximum();
   h_draw_st_fit->SetMaximum(max*1.1);
   
   legf->Draw();
   
   std::string fsave = "RESULT.eps";
   c1->Print(fsave.c_str());
   delete legf;
   delete h_draw_st_fit;
   delete c1;
}

// uncertainty band creation
TGraphAsymmErrors* CFIT::cfit::makeErrorBand(TH1D* tot, TH1D* plus, TH1D* minus)
{   
   int nbins = tot->GetNbinsX();
   
   const int s = 1000;
   
   double x[s];
   double xerr[s];
   double y[s];
   double yp[s];
   double ym[s];
   
   for(int bin =1; bin<=nbins; ++bin)
     {    
	int index = bin-1;
	xerr[index] = tot->GetBinWidth(bin)/2.0;
	x[index] = tot->GetBinLowEdge(bin) + xerr[index];
	
	y[index] = tot->GetBinContent(bin);
	yp[index] = plus->GetBinContent(bin)-y[index];
	ym[index] = y[index]-minus->GetBinContent(bin);
	if(y[index] - ym[index] < 0) ym[index] = y[index];
     }
   
   TGraphAsymmErrors* error = new TGraphAsymmErrors(nbins, x, y, xerr, xerr, ym, yp);
   
   return error;
}

void CFIT::cfit::removeMCSys(TH1D *hnom,TH1D *hsysUp,TH1D *hsysLow)
{
   int nb = hnom->GetXaxis()->GetNbins();
   for(int ib=1;ib<=nb;ib++)
     {
        double vb0 = hnom->GetBinContent(ib);
        double ve0 = hnom->GetBinError(ib);
        double vbUp = hsysUp->GetBinContent(ib);

        if( vbUp-vb0 > vb0 )
          {
             hsysUp->SetBinContent(ib,2*vb0);
          }
     }
}

void CFIT::cfit::silentMode() { verb = 0; }

double CFIT::cfit::getCHISQ() { return CHISQ; }

double CFIT::cfit::getPAR(int idx) 
{ 
   if( idx >= nT )
     {
	std::cout << "The number of templates used is " << nT << std::endl;
	exit(1);
     }   
   
   if( VALID ) return PAR[idx];
   else return -1.;
}

double CFIT::cfit::getERR(int idx)
{ 
   if( idx >= nT )
     {
	std::cout << "The number of templates used is " << nT << std::endl;
	exit(1);
     }   
   
   if( VALID ) return ERR[idx];
   else return -1.;
}
