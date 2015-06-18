#include "../cfit.h"

void getResults(CFIT::cfit *cf,float *par,float *err);

void test()
{
   gROOT->SetBatch();
   
   gSystem->Load("../libCFIT.so");

   float par[100];
   float err[100];
   
   CFIT::cfit *cf = new CFIT::cfit("Fit variable");
   cf->SetOptimization(OPT_MORPH_SGN_SIGMA);
   cf->SetMorphing(OPTMORPH_GEOMETRIC);
   cf->ProducePlots(1);
   
   cf->SetInputFile("test.root");

   cf->AddSys("SYS1","_sys1_down","_sys1_up");
   cf->AddSys("SYS2","_sys2_down","_sys2_up");

   cf->SetMatrixOption("WRITE");
   
   cf->SetData("h_data");
   
   cf->AddTemplate("template1","h_mc1");
   cf->AddTemplate("template2","h_mc2");
   cf->AddTemplate("template3","h_mc3");

   cf->Run();
   
   getResults(cf,par,err);
   float chi2 = cf->GetChisq();
   
   std::cout << "chi2=" << chi2 << std::endl;
   for(int i=0;i<3;i++)
     {
	std::cout << par[i] << " +- " << err[i] << std::endl;
     }   

   delete cf;
   
   gApplication->Terminate();
}

void getResults(CFIT::cfit *cf,float *par,float *err)
{   
   float chi2 = cf->GetChisq();
   int nPar = cf->GetNPar();
   for(int i=0;i<nPar;i++)
     {
	par[i] = cf->GetPar(i);
	err[i] = cf->GetParErr(i);	
     }
}
