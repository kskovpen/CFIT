import sys
sys.path.insert(0, '../')
from cfit import *


def testSF():
   par = []
   err = []
   par_tag = []
   err_tag = []

   cf = cfit("Fit variable")
   cf.SetOptimization(OPT_MORPH_SGN_SIGMA)
   cf.SetMorphing(OPTMORPH_GEOMETRIC)
   cf.ProducePlots(True)
   
   cf.SetInputFile("test.root")

   cf.AddSys("SYS1","_sys1_down","_sys1_up")
   cf.AddSys("SYS2","_sys2_down","_sys2_up")

   cf.SetMatrixOption("WRITE")
   
   cf.SetData("h_data")
   cf.SetDataTag("h_data_tag")
   cf.SetDataUntag("h_data_untag")
   
   cf.AddTemplate("template1","h_mc1",2)
   cf.AddTemplate("template2","h_mc2",3)
   cf.AddTemplate("template3","h_mc3",4)

   cf.AddTemplateTag("template1","h_mc1_tag",2)
   cf.AddTemplateTag("template2","h_mc2_tag",3)
   cf.AddTemplateTag("template3","h_mc3_tag",4)

   cf.AddTemplateUntag("template1","h_mc1_untag",2)
   cf.AddTemplateUntag("template2","h_mc2_untag",3)
   cf.AddTemplateUntag("template3","h_mc3_untag",4)
   
   cf.Run()   
   for i in range(cf.GetNPar()):
      par.append(cf.GetPar(i))
      err.append(cf.GetParErr(i))
   
   chi2 = cf.GetChisq()
   ndata = cf.GetNData()
   nmc1 = cf.GetNTemplate("template1")
   nmc = cf.GetNTemplate("template1")+cf.GetNTemplate("template2")+cf.GetNTemplate("template3")

   cf.Run("tag")   
   for i in range(cf.GetNPar()):
      par_tag.append(cf.GetPar(i))
      err_tag.append(cf.GetParErr(i))
   chi2_tag = cf.GetChisq();
   ndata_tag = cf.GetNData();
   nmc1_tag = cf.GetNTemplate("template1");
   nmc_tag = cf.GetNTemplate("template1")+cf.GetNTemplate("template2")+cf.GetNTemplate("template3");

   fr = nmc1/nmc;
   fr_tag = nmc1_tag/nmc_tag;
   
   effMC = nmc1_tag/nmc1;
   effDATA = nmc_tag/nmc*par_tag[0]/par[0]*fr_tag/fr;
   
   print "effMC = ",  effMC
   print "effDATA = ", effDATA 
   print "sf = ", effDATA/effMC

   
   print par
   print err
   par, err = [], []
   # perform statistical variation
   cf.SetMatrixOption("READ");
   cf.SetStatVariation(667);
   cf.ProducePlots(False);   
   cf.Run();
   for i in range(cf.GetNPar()):
      par.append(cf.GetPar(i))
      err.append(cf.GetParErr(i))
   cf.SetStatVariation(667);
   cf.Run("tag");
   for i in range(cf.GetNPar()):
      par.append(cf.GetPar(i))
      err.append(cf.GetParErr(i))
   
   print par
   print err
   par, err = [], []
   # do the calculation of SF here again ....
   # perform systematic variation
   cf.SetMatrixOption("READ");
   cf.SetSysVariation("_sys1_up");
   cf.Run();
   for i in range(cf.GetNPar()):
      par.append(cf.GetPar(i))
      err.append(cf.GetParErr(i))
   cf.SetSysVariation("_sys1_up");
   cf.Run("tag");
   for i in range(cf.GetNPar()):
      par.append(cf.GetPar(i))
      err.append(cf.GetParErr(i))

   
   print par
   print err

testSF()