import sys
sys.path.insert(0, '../')
from cfit import *

def test():
  cf = cfit("Fit variable")
  cf.SetOptimization(OPT_MORPH_SGN_SIGMA)
  cf.SetMorphing(OPTMORPH_GEOMETRIC)
  cf.ProducePlots(True)
   
  cf.SetInputFile("test.root")

  cf.AddSys("SYS1","_sys1_down","_sys1_up")
  cf.AddSys("SYS2","_sys2_down","_sys2_up")

  cf.SetMatrixOption("WRITE")
   
  cf.SetData("h_data")
   
  cf.AddTemplate("template1","h_mc1",2)
  cf.AddTemplate("template2","h_mc2",3)
  cf.AddTemplate("template3","h_mc3",4)

  cf.Run()
   
  chi2 = cf.GetChisq();   
  print "chi2 =", chi2
  
  for i in range(cf.GetNPar()):
    par_i = cf.GetPar(i)
    err_i = cf.GetParErr(i)
    print par_i, "+-", err_i

test()
