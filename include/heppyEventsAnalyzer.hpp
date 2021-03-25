//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 11 10:41:01 2020 by ROOT version 6.16/00
// from TTree events/NtupleProducer
// found on file: MC_tW_top.root
//////////////////////////////////////////////////////////


#ifndef heppyEventsAnalyzer_h
#define heppyEventsAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <string>
#include <TClonesArray.h>
#include <TObject.h>
using namespace std;

// Header file for the classes stored in the TTree if any.

class heppyEventsAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   const char *	 InputFile;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           is_data;
   Int_t           NumberOfLeptons;
   Int_t           natureLepton;
   Int_t           run;
   Int_t           region;
   ULong64_t       event;
   Double_t        rho;
   Int_t           lumi;
   Int_t           n_pu;
   Int_t           n_up;
   Double_t        weight_generator;
   Int_t           n_pv;
   Double_t        weight_pu;
   Double_t        weight_sfm_id;
   Double_t        weight;
   Double_t        weight_sfm_trig_isomu27;
   Double_t        weight_sfb;
   Double_t        weight_sfe_reco;
   Double_t        weight_sfm_trig_mu50;
   Double_t        weight_sfm_iso;
   Double_t        weight_sf_em_trig;
   Double_t        weight_sfe_id;
   Double_t        j3_pumva;
   Double_t        j1_pt;
   Int_t           n_jets_pt40;
   Double_t        j2_rawf;
   Double_t        j3_energy;
   Double_t        j3_rawf;
   Double_t        j1_rawf;
   Double_t        j3_flavour_parton;
   Double_t        j2_flavour_hadron;
   Double_t        j3_phi;
   Double_t        j3_eta;
   Double_t        j3_mass;
   Double_t        j1_eta;
   Double_t        j1_phi;
   Double_t        j2_energy;
   Double_t        j1_energy;
   Double_t        j2_flavour_parton;
   Double_t        j1_flavour_hadron;
   Double_t        j3_flavour_hadron;
   Double_t        j2_phi;
   Double_t        j1_flavour_parton;
   Double_t        j2_mass;
   Double_t        j2_pt;
   Double_t        j1_pumva;
   Double_t        j2_pumva;
   Double_t        j3_pt;
   Double_t        j1_mass;
   Double_t        j2_eta;
   Double_t        b2_phi;
   Double_t        b1_pt;
   Double_t        b3_flavour_hadron;
   Double_t        b1_flavour_parton;
   Double_t        b3_eta;
   Double_t        b3_mass;
   Double_t        b3_pumva;
   Double_t        b1_energy;
   Double_t        b3_energy;
   Double_t        b2_mass;
   Double_t        b2_eta;
   Double_t        b2_flavour_hadron;
   Double_t        b3_pt;
   Double_t        b1_pumva;
   Double_t        b3_flavour_parton;
   Double_t        b2_rawf;
   Double_t        b1_phi;
   Double_t        b2_pt;
   Double_t        b3_phi;
   Double_t        b1_flavour_hadron;
   Double_t        b2_pumva;
   Int_t           n_bjets;
   Double_t        b1_mass;
   Double_t        b3_rawf;
   Double_t        b1_eta;
   Double_t        b1_rawf;
   Double_t        b2_energy;
   Double_t        b2_flavour_parton;
   Double_t        eta_elec;
   Double_t        phi_elec;
   Double_t        energy_elec;
   Double_t        pt_elec;
   Double_t        m_elec;
   Double_t        iso_elec;
   Double_t        q_elec;
   Double_t        energy_muon;
   Double_t        pt_muon;
   Double_t        eta_muon;
   Double_t        q_muon;
   Double_t        phi_muon;
   Double_t        m_muon;
   Double_t        iso_muon;
   Double_t        metphi;
   Double_t        met;
   Double_t        trg_double_muon_mu17m3_fired;
   Double_t        trg_electron_ele35_fired;
   Double_t        trg_electron_ele40_fired;
   Double_t        trg_muon_mu24tk_fired;
   Double_t        trg_electron_ele32doubleEG_fired;
   Double_t        trg_double_muon_mu17_fired;
   Double_t        trg_muon_mu27_fired;
   Double_t        trg_electron_ele32_fired;
   Double_t        trg_muon_mu24eta21_fired;
   Double_t        trg_muon_mu24_fired;
   Double_t        trg_electron_ele38_fired;

   // List of branches
   TBranch        *b_is_data;   //!
   TBranch        *b_NumberOfLeptons;   //!
   TBranch        *b_natureLepton;   //!
   TBranch        *b_run;   //!
   TBranch        *b_region;   //!
   TBranch        *b_event;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_n_pu;   //!
   TBranch        *b_n_up;   //!
   TBranch        *b_weight_generator;   //!
   TBranch        *b_n_pv;   //!
   TBranch        *b_weight_pu;   //!
   TBranch        *b_weight_sfm_id;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_weight_sfm_trig_isomu27;   //!
   TBranch        *b_weight_sfb;   //!
   TBranch        *b_weight_sfe_reco;   //!
   TBranch        *b_weight_sfm_trig_mu50;   //!
   TBranch        *b_weight_sfm_iso;   //!
   TBranch        *b_weight_sf_em_trig;   //!
   TBranch        *b_weight_sfe_id;   //!
   TBranch        *b_j3_pumva;   //!
   TBranch        *b_j1_pt;   //!
   TBranch        *b_n_jets_pt40;   //!
   TBranch        *b_j2_rawf;   //!
   TBranch        *b_j3_energy;   //!
   TBranch        *b_j3_rawf;   //!
   TBranch        *b_j1_rawf;   //!
   TBranch        *b_j3_flavour_parton;   //!
   TBranch        *b_j2_flavour_hadron;   //!
   TBranch        *b_j3_phi;   //!
   TBranch        *b_j3_eta;   //!
   TBranch        *b_j3_mass;   //!
   TBranch        *b_j1_eta;   //!
   TBranch        *b_j1_phi;   //!
   TBranch        *b_j2_energy;   //!
   TBranch        *b_j1_energy;   //!
   TBranch        *b_j2_flavour_parton;   //!
   TBranch        *b_j1_flavour_hadron;   //!
   TBranch        *b_j3_flavour_hadron;   //!
   TBranch        *b_j2_phi;   //!
   TBranch        *b_j1_flavour_parton;   //!
   TBranch        *b_j2_mass;   //!
   TBranch        *b_j2_pt;   //!
   TBranch        *b_j1_pumva;   //!
   TBranch        *b_j2_pumva;   //!
   TBranch        *b_j3_pt;   //!
   TBranch        *b_j1_mass;   //!
   TBranch        *b_j2_eta;   //!
   TBranch        *b_b2_phi;   //!
   TBranch        *b_b1_pt;   //!
   TBranch        *b_b3_flavour_hadron;   //!
   TBranch        *b_b1_flavour_parton;   //!
   TBranch        *b_b3_eta;   //!
   TBranch        *b_b3_mass;   //!
   TBranch        *b_b3_pumva;   //!
   TBranch        *b_b1_energy;   //!
   TBranch        *b_b3_energy;   //!
   TBranch        *b_b2_mass;   //!
   TBranch        *b_b2_eta;   //!
   TBranch        *b_b2_flavour_hadron;   //!
   TBranch        *b_b3_pt;   //!
   TBranch        *b_b1_pumva;   //!
   TBranch        *b_b3_flavour_parton;   //!
   TBranch        *b_b2_rawf;   //!
   TBranch        *b_b1_phi;   //!
   TBranch        *b_b2_pt;   //!
   TBranch        *b_b3_phi;   //!
   TBranch        *b_b1_flavour_hadron;   //!
   TBranch        *b_b2_pumva;   //!
   TBranch        *b_n_bjets;   //!
   TBranch        *b_b1_mass;   //!
   TBranch        *b_b3_rawf;   //!
   TBranch        *b_b1_eta;   //!
   TBranch        *b_b1_rawf;   //!
   TBranch        *b_b2_energy;   //!
   TBranch        *b_b2_flavour_parton;   //!
   TBranch        *b_eta_elec;   //!
   TBranch        *b_phi_elec;   //!
   TBranch        *b_energy_elec;   //!
   TBranch        *b_pt_elec;   //!
   TBranch        *b_m_elec;   //!
   TBranch        *b_iso_elec;   //!
   TBranch        *b_q_elec;   //!
   TBranch        *b_energy_muon;   //!
   TBranch        *b_pt_muon;   //!
   TBranch        *b_eta_muon;   //!
   TBranch        *b_q_muon;   //!
   TBranch        *b_phi_muon;   //!
   TBranch        *b_m_muon;   //!
   TBranch        *b_iso_muon;   //!
   TBranch        *b_metphi;   //!
   TBranch        *b_met;   //!
   TBranch        *b_trg_electron_ele35_fired;   //!
   TBranch        *b_trg_electron_ele40_fired;   //!
   TBranch        *b_trg_muon_electron_mu12ele23DZ_fired;   //!
   TBranch        *b_trg_muon_mu24tk_fired;   //!
   TBranch        *b_trg_electron_ele32doubleEG_fired;   //!
   TBranch        *b_trg_double_muon_mu17_fired;   //!
   TBranch        *b_trg_muon_mu27_fired;   //!
   TBranch        *b_trg_electron_ele32_fired;   //!
   TBranch        *b_trg_muon_mu24eta21_fired;   //!
   TBranch        *b_trg_muon_mu24_fired;   //!
   TBranch        *b_trg_electron_ele38_fired;   //!

   heppyEventsAnalyzer(TTree *tree=0);
   virtual ~heppyEventsAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef heppyEventsAnalyzer_cxx
heppyEventsAnalyzer::heppyEventsAnalyzer(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(InputFile);
      if (!f || !f->IsOpen()) {
         f = new TFile(InputFile);
      }
      f->GetObject("events",tree);
      cout<<InputFile<<endl;
   }
   Init(tree);
}

heppyEventsAnalyzer::~heppyEventsAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t heppyEventsAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t heppyEventsAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}
void heppyEventsAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("is_data", &is_data, &b_is_data);
   fChain->SetBranchAddress("NumberOfLeptons", &NumberOfLeptons, &b_NumberOfLeptons);
   fChain->SetBranchAddress("natureLepton", &natureLepton, &b_natureLepton);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("region", &region, &b_region);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("n_pu", &n_pu, &b_n_pu);
   fChain->SetBranchAddress("n_up", &n_up, &b_n_up);
   fChain->SetBranchAddress("weight_generator", &weight_generator, &b_weight_generator);
   fChain->SetBranchAddress("n_pv", &n_pv, &b_n_pv);
   fChain->SetBranchAddress("weight_pu", &weight_pu, &b_weight_pu);
   fChain->SetBranchAddress("weight_sfm_id", &weight_sfm_id, &b_weight_sfm_id);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("weight_sfm_trig_isomu27", &weight_sfm_trig_isomu27, &b_weight_sfm_trig_isomu27);
   fChain->SetBranchAddress("weight_sfb", &weight_sfb, &b_weight_sfb);
   fChain->SetBranchAddress("weight_sfe_reco", &weight_sfe_reco, &b_weight_sfe_reco);
   fChain->SetBranchAddress("weight_sfm_trig_mu50", &weight_sfm_trig_mu50, &b_weight_sfm_trig_mu50);
   fChain->SetBranchAddress("weight_sfm_iso", &weight_sfm_iso, &b_weight_sfm_iso);
   fChain->SetBranchAddress("weight_sf_em_trig", &weight_sf_em_trig, &b_weight_sf_em_trig);
   fChain->SetBranchAddress("weight_sfe_id", &weight_sfe_id, &b_weight_sfe_id);
   fChain->SetBranchAddress("j3_pumva", &j3_pumva, &b_j3_pumva);
   fChain->SetBranchAddress("j1_pt", &j1_pt, &b_j1_pt);
   fChain->SetBranchAddress("n_jets_pt40", &n_jets_pt40, &b_n_jets_pt40);
   fChain->SetBranchAddress("j2_rawf", &j2_rawf, &b_j2_rawf);
   fChain->SetBranchAddress("j3_energy", &j3_energy, &b_j3_energy);
   fChain->SetBranchAddress("j3_rawf", &j3_rawf, &b_j3_rawf);
   fChain->SetBranchAddress("j1_rawf", &j1_rawf, &b_j1_rawf);
   fChain->SetBranchAddress("j3_flavour_parton", &j3_flavour_parton, &b_j3_flavour_parton);
   fChain->SetBranchAddress("j2_flavour_hadron", &j2_flavour_hadron, &b_j2_flavour_hadron);
   fChain->SetBranchAddress("j3_phi", &j3_phi, &b_j3_phi);
   fChain->SetBranchAddress("j3_eta", &j3_eta, &b_j3_eta);
   fChain->SetBranchAddress("j3_mass", &j3_mass, &b_j3_mass);
   fChain->SetBranchAddress("j1_eta", &j1_eta, &b_j1_eta);
   fChain->SetBranchAddress("j1_phi", &j1_phi, &b_j1_phi);
   fChain->SetBranchAddress("j2_energy", &j2_energy, &b_j2_energy);
   fChain->SetBranchAddress("j1_energy", &j1_energy, &b_j1_energy);
   fChain->SetBranchAddress("j2_flavour_parton", &j2_flavour_parton, &b_j2_flavour_parton);
   fChain->SetBranchAddress("j1_flavour_hadron", &j1_flavour_hadron, &b_j1_flavour_hadron);
   fChain->SetBranchAddress("j3_flavour_hadron", &j3_flavour_hadron, &b_j3_flavour_hadron);
   fChain->SetBranchAddress("j2_phi", &j2_phi, &b_j2_phi);
   fChain->SetBranchAddress("j1_flavour_parton", &j1_flavour_parton, &b_j1_flavour_parton);
   fChain->SetBranchAddress("j2_mass", &j2_mass, &b_j2_mass);
   fChain->SetBranchAddress("j2_pt", &j2_pt, &b_j2_pt);
   fChain->SetBranchAddress("j1_pumva", &j1_pumva, &b_j1_pumva);
   fChain->SetBranchAddress("j2_pumva", &j2_pumva, &b_j2_pumva);
   fChain->SetBranchAddress("j3_pt", &j3_pt, &b_j3_pt);
   fChain->SetBranchAddress("j1_mass", &j1_mass, &b_j1_mass);
   fChain->SetBranchAddress("j2_eta", &j2_eta, &b_j2_eta);
   fChain->SetBranchAddress("b2_phi", &b2_phi, &b_b2_phi);
   fChain->SetBranchAddress("b1_pt", &b1_pt, &b_b1_pt);
   fChain->SetBranchAddress("b3_flavour_hadron", &b3_flavour_hadron, &b_b3_flavour_hadron);
   fChain->SetBranchAddress("b1_flavour_parton", &b1_flavour_parton, &b_b1_flavour_parton);
   fChain->SetBranchAddress("b3_eta", &b3_eta, &b_b3_eta);
   fChain->SetBranchAddress("b3_mass", &b3_mass, &b_b3_mass);
   fChain->SetBranchAddress("b3_pumva", &b3_pumva, &b_b3_pumva);
   fChain->SetBranchAddress("b1_energy", &b1_energy, &b_b1_energy);
   fChain->SetBranchAddress("b3_energy", &b3_energy, &b_b3_energy);
   fChain->SetBranchAddress("b2_mass", &b2_mass, &b_b2_mass);
   fChain->SetBranchAddress("b2_eta", &b2_eta, &b_b2_eta);
   fChain->SetBranchAddress("b2_flavour_hadron", &b2_flavour_hadron, &b_b2_flavour_hadron);
   fChain->SetBranchAddress("b3_pt", &b3_pt, &b_b3_pt);
   fChain->SetBranchAddress("b1_pumva", &b1_pumva, &b_b1_pumva);
   fChain->SetBranchAddress("b3_flavour_parton", &b3_flavour_parton, &b_b3_flavour_parton);
   fChain->SetBranchAddress("b2_rawf", &b2_rawf, &b_b2_rawf);
   fChain->SetBranchAddress("b1_phi", &b1_phi, &b_b1_phi);
   fChain->SetBranchAddress("b2_pt", &b2_pt, &b_b2_pt);
   fChain->SetBranchAddress("b3_phi", &b3_phi, &b_b3_phi);
   fChain->SetBranchAddress("b1_flavour_hadron", &b1_flavour_hadron, &b_b1_flavour_hadron);
   fChain->SetBranchAddress("b2_pumva", &b2_pumva, &b_b2_pumva);
   fChain->SetBranchAddress("n_bjets", &n_bjets, &b_n_bjets);
   fChain->SetBranchAddress("b1_mass", &b1_mass, &b_b1_mass);
   fChain->SetBranchAddress("b3_rawf", &b3_rawf, &b_b3_rawf);
   fChain->SetBranchAddress("b1_eta", &b1_eta, &b_b1_eta);
   fChain->SetBranchAddress("b1_rawf", &b1_rawf, &b_b1_rawf);
   fChain->SetBranchAddress("b2_energy", &b2_energy, &b_b2_energy);
   fChain->SetBranchAddress("b2_flavour_parton", &b2_flavour_parton, &b_b2_flavour_parton);
   fChain->SetBranchAddress("eta_elec", &eta_elec, &b_eta_elec);
   fChain->SetBranchAddress("phi_elec", &phi_elec, &b_phi_elec);
   fChain->SetBranchAddress("energy_elec", &energy_elec, &b_energy_elec);
   fChain->SetBranchAddress("pt_elec", &pt_elec, &b_pt_elec);
   fChain->SetBranchAddress("m_elec", &m_elec, &b_m_elec);
   fChain->SetBranchAddress("iso_elec", &iso_elec, &b_iso_elec);
   fChain->SetBranchAddress("q_elec", &q_elec, &b_q_elec);
   fChain->SetBranchAddress("energy_muon", &energy_muon, &b_energy_muon);
   fChain->SetBranchAddress("pt_muon", &pt_muon, &b_pt_muon);
   fChain->SetBranchAddress("eta_muon", &eta_muon, &b_eta_muon);
   fChain->SetBranchAddress("q_muon", &q_muon, &b_q_muon);
   fChain->SetBranchAddress("phi_muon", &phi_muon, &b_phi_muon);
   fChain->SetBranchAddress("m_muon", &m_muon, &b_m_muon);
   fChain->SetBranchAddress("iso_muon", &iso_muon, &b_iso_muon);
   fChain->SetBranchAddress("metphi", &metphi, &b_metphi);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("trg_electron_ele35_fired", &trg_electron_ele35_fired, &b_trg_electron_ele35_fired);
   fChain->SetBranchAddress("trg_electron_ele40_fired", &trg_electron_ele40_fired, &b_trg_electron_ele40_fired);
   fChain->SetBranchAddress("trg_muon_mu24tk_fired", &trg_muon_mu24tk_fired, &b_trg_muon_mu24tk_fired);
   fChain->SetBranchAddress("trg_electron_ele32doubleEG_fired", &trg_electron_ele32doubleEG_fired, &b_trg_electron_ele32doubleEG_fired);
   fChain->SetBranchAddress("trg_muon_mu27_fired", &trg_muon_mu27_fired, &b_trg_muon_mu27_fired);
   fChain->SetBranchAddress("trg_electron_ele32_fired", &trg_electron_ele32_fired, &b_trg_electron_ele32_fired);
   fChain->SetBranchAddress("trg_muon_mu24eta21_fired", &trg_muon_mu24eta21_fired, &b_trg_muon_mu24eta21_fired);
   fChain->SetBranchAddress("trg_muon_mu24_fired", &trg_muon_mu24_fired, &b_trg_muon_mu24_fired);
   fChain->SetBranchAddress("trg_electron_ele38_fired", &trg_electron_ele38_fired, &b_trg_electron_ele38_fired);
   Notify();
}

Bool_t heppyEventsAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void heppyEventsAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t heppyEventsAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analyzer_cxx
