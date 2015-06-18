#include "../cfit.h"

void getResults(CFIT::cfit *cf,float *par,float *err);

void testSF()
{
   gROOT->SetBatch();
   
   gSystem->Load("../libCFIT.so");

   float par[100];
   float err[100];
   float par_tag[100];
   float err_tag[100];
   
   CFIT::cfit *cf = new CFIT::cfit("Fit variable");
   cf->SetOptimization(OPT_MORPH_SGN_SIGMA);
   cf->SetMorphing(OPTMORPH_GEOMETRIC);
   cf->ProducePlots(1);
   
   cf->SetInputFile("test.root");

   cf->AddSys("SYS1","_sys1_down","_sys1_up");
   cf->AddSys("SYS2","_sys2_down","_sys2_up");

   cf->SetMatrixOption("WRITE");
   
   cf->SetData("h_data");
   cf->SetDataTag("h_data_tag");
   cf->SetDataUntag("h_data_untag");
   
   cf->AddTemplate("template1","h_mc1");
   cf->AddTemplate("template2","h_mc2");
   cf->AddTemplate("template3","h_mc3");

   cf->AddTemplateTag("template1","h_mc1_tag");
   cf->AddTemplateTag("template2","h_mc2_tag");
   cf->AddTemplateTag("template3","h_mc3_tag");

   cf->AddTemplateUntag("template1","h_mc1_untag");
   cf->AddTemplateUntag("template2","h_mc2_untag");
   cf->AddTemplateUntag("template3","h_mc3_untag");
   
   cf->Run();   
   getResults(cf,par,err);
   float chi2 = cf->GetChisq();
   float ndata = cf->GetNData();
   float nmc1 = cf->GetNTemplate("template1");
   float nmc = cf->GetNTemplate("template1")+cf->GetNTemplate("template2")+cf->GetNTemplate("template3");

   cf->Run("tag");   
   getResults(cf,par_tag,err_tag);
   float chi2_tag = cf->GetChisq();
   float ndata_tag = cf->GetNData();
   float nmc1_tag = cf->GetNTemplate("template1");
   float nmc_tag = cf->GetNTemplate("template1")+cf->GetNTemplate("template2")+cf->GetNTemplate("template3");

   float fr = nmc1/nmc;
   float fr_tag = nmc1_tag/nmc_tag;
   
   float effMC = nmc1_tag/nmc1;
   float effDATA = nmc_tag/nmc*par_tag[0]/par[0]*fr_tag/fr;
   
   std::cout << "effMC = " << effMC << std::endl;
   std::cout << "effDATA = " << effDATA << std::endl;
   std::cout << "sf = " << effDATA/effMC << std::endl;

   // perform statistical variation
   cf->SetMatrixOption("READ");
   cf->SetStatVariation(667);
   cf->ProducePlots(0);
   
   cf->Run();
   getResults(cf,par,err);

   cf->Run("tag");
   getResults(cf,par,err);
   
   // do the calculation of SF here again ....
   
   // perform systematic variation
   cf->SetMatrixOption("READ");
   cf->SetSysVariation("_sys1_up");

   cf->Run();
   getResults(cf,par,err);

   cf->Run("tag");
   getResults(cf,par,err);
   
   // do the calculation of SF here again ....
   
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
