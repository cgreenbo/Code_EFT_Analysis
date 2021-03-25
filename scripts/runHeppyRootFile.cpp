
#include "../include/heppyEventsAnalyzer.hpp"
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <sstream>
#include <iostream>

using namespace std;


int main()
{

  int nSamples =2;
  string* name = new string[nSamples];
  /*name[0] = "MC_t_top.root";
  name[1] = "MC_t_antitop.root";
  name[2] = "MC_DY_50.root";
  name[3] = "MC_DY_1050.root";
  name[4] = "MC_semilep.root";
  name[5] = "MC_tW_antitop.root";
  name[6] = "MC_tW_top.root";
  name[7] = "MC_W+Jets.root";
  name[8] = "MC_dilep.root";
  name[9] = "DATA_region.root";
  name[10] = "multijet_elec_new_trigger.root";
  name[11] = "multijet_muon_retry.root";
  name[12] = "MC_QCD_EM_20to30.root";
  name[13] = "MC_QCD_EM_30to50.root";
  name[14] = "MC_QCD_EM_50to80.root";
  name[15] = "MC_QCD_EM_80to120.root";
  name[16] = "MC_QCD_EM_120to170.root";
  name[17] = "MC_QCD_EM_170to300.root";
  name[18] = "MC_QCD_EM_300toInf.root";
  name[19] = "MC_QCD_Mu_15to20.root";
  name[20] = "MC_QCD_Mu_20to30.root";
  name[21] = "MC_QCD_Mu_30to50.root";
  name[22] = "MC_QCD_Mu_50to80.root";
  name[23] = "MC_QCD_Mu_80to120.root";
  name[24] = "MC_QCD_Mu_120to170.root";
  name[25] = "MC_QCD_Mu_170to300.root";
  name[26] = "MC_QCD_Mu_300to470.root";
  name[27] = "MC_QCD_Mu_470to600.root";
  name[28] = "MC_QCD_Mu_600to800.root";
  name[29] = "MC_QCD_Mu_800to1000.root";
  name[30] = "MC_QCD_Mu_1000toInf.root";*/

  name[0] = "multijet_elec_0_muon_iso.root";
  name[1] = "multijet_muon_0_elec_iso.root";
  /*name[0] = "MC_t_top_QCD.root";
  name[1] = "MC_t_antitop_QCD.root";
  name[2] = "MC_DY_50_QCD.root";
  name[3] = "MC_DY_1050_QCD.root";
  name[4] = "MC_semilep_QCD.root";
  name[5] = "MC_tW_antitop_QCD.root";
  name[6] = "MC_tW_top_QCD.root";
  name[7] = "MC_W+Jets_QCD.root";
  name[8] = "MC_dilep_QCD.root";
  name[9] = "DATA_QCD.root";
  name[10] = "MC_QCD_EM_20to30_QCD.root";
  name[11] = "MC_QCD_EM_30to50_QCD.root";
  name[12] = "MC_QCD_EM_50to80_QCD.root";
  name[13] = "MC_QCD_EM_80to120_QCD.root";
  name[14] = "MC_QCD_EM_120to170_QCD.root";
  name[15] = "MC_QCD_EM_170to300_QCD.root";
  name[16] = "MC_QCD_EM_300toInf_QCD.root";
  name[17] = "MC_QCD_Mu_15to20_QCD.root";
  name[18] = "MC_QCD_Mu_20to30_QCD.root";
  name[19] = "MC_QCD_Mu_30to50_QCD.root";
  name[20] = "MC_QCD_Mu_50to80_QCD.root";
  name[21] = "MC_QCD_Mu_80to120_QCD.root";
  name[22] = "MC_QCD_Mu_120to170_QCD.root";
  name[23] = "MC_QCD_Mu_170to300_QCD.root";
  name[24] = "MC_QCD_Mu_300to470_QCD.root";
  name[25] = "MC_QCD_Mu_470to600_QCD.root";
  name[26] = "MC_QCD_Mu_600to800_QCD.root";
  name[27] = "MC_QCD_Mu_800to1000_QCD.root";
  name[28] = "MC_QCD_Mu_1000toInf_QCD.root";*/





  TFile** fInput = new TFile*[nSamples];
  TTree** tInput = new TTree*[nSamples];

  heppyEventsAnalyzer** HeppyEventsAnalyzer = new heppyEventsAnalyzer*[nSamples];

    for (int i = 0; i<nSamples; i++)
    {
      cout<<"Treating rootfile "+name[i]<<endl;
      string inputName = "data/heppy/"+name[i];
      string outputName = "output_"+name[i];

      fInput[i] = new TFile(inputName.c_str(),"READ");
      tInput[i] = (TTree*) fInput[i]->Get("events");
      HeppyEventsAnalyzer[i] = new heppyEventsAnalyzer(tInput[i]);
      HeppyEventsAnalyzer[i]->Loop();
      string commandline = "mv output.root data/heppy/output/" + outputName;
      system(commandline.c_str());
    }


return 0;
}
