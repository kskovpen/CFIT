void test()
{
   gROOT->SetBatch();
   
   // load CFIT library
   gSystem->Load("libCFIT.so");

   // create an instance of CFIT
   CFIT::cfit cf("input.txt");
//   cf.silentMode();
   
   // run CFIT
   cf.run();
 
   // user methods
//   std::cout << cf.getCHISQ() << std::endl;
//   for(int i=0;i<3;i++)
//     std::cout << cf.getPAR(i) << " +- " << cf.getERR(i) << std::endl;
   
   gApplication->Terminate();
}
