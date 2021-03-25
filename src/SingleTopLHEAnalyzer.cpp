#define SingleTopLHEAnalyzer_cxx
#include "../include/SingleTopLHEAnalyzer.hpp"

#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

void SingleTopLHEAnalyzer::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L SingleTopLHEAnalyzer.C
//      Root > SingleTopLHEAnalyzer t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   float weight;
   float sinTheta, cosTheta;
   float sinThetaStar, cosThetaStar;
   float sinPhiStar, cosPhiStar, PhiStar;
   float lepton_E_Wframe;
   float top_pt, W_pt, lepton_pt;
   float top_mass, W_mass, W_transverse_mass;
   float nature_lepton;
   float weight_SM;
   float weight_ctwi_m5, weight_ctwi_m2, weight_ctwi_m1, weight_ctwi_p1, weight_ctwi_p2, weight_ctwi_p5;
   float weight_cbwi_m5, weight_cbwi_m2, weight_cbwi_m1, weight_cbwi_p1, weight_cbwi_p2, weight_cbwi_p5;
   
   TFile* fOutput = new TFile("output.root","RECREATE");
   TTree* tOutput = new TTree("Tree","Tree");

   tOutput->Branch("nature_lepton",&nature_lepton,"nature_lepton/F");
   tOutput->Branch("weight",&weight,"weight/F");
   tOutput->Branch("cosTheta",&cosTheta,"cosTheta/F");
   tOutput->Branch("sinTheta",&sinTheta,"sinTheta/F");
   tOutput->Branch("cosThetaStar",&cosThetaStar,"cosThetaStar/F");
   tOutput->Branch("sinThetaStar",&sinThetaStar,"sinThetaStar/F");
   tOutput->Branch("cosPhiStar",&cosPhiStar,"cosPhiStar/F");
   tOutput->Branch("sinPhiStar",&sinPhiStar,"sinPhiStar/F");
   tOutput->Branch("PhiStar",&PhiStar,"PhiStar/F");
   tOutput->Branch("lepton_E_Wframe",&lepton_E_Wframe,"lepton_E_Wframe/F");
   tOutput->Branch("top_pt",&top_pt,"top_pt/F");
   tOutput->Branch("W_pt",&W_pt,"W_pt/F");
   tOutput->Branch("lepton_pt",&lepton_pt,"lepton_pt/F");
   tOutput->Branch("W_mass", &W_mass, "W_mass/F");
   tOutput->Branch("top_mass", &top_mass, "top_mass/F");
   tOutput->Branch("W_transverse_mass", &W_transverse_mass, "W_transverse_mass/F");

   tOutput->Branch("weight_SM",&weight_SM,"weight_SM/F");

   //tOutput->Branch("weight_ctwi_m5",&weight_ctwi_m5,"weight_ctwi_m5/F");
   tOutput->Branch("weight_ctwi_m2",&weight_ctwi_m2,"weight_ctwi_m2/F");
   tOutput->Branch("weight_ctwi_m1",&weight_ctwi_m1,"weight_ctwi_m1/F");
   tOutput->Branch("weight_ctwi_p1",&weight_ctwi_p1,"weight_ctwi_p1/F");
   tOutput->Branch("weight_ctwi_p2",&weight_ctwi_p2,"weight_ctwi_p2/F");
   //tOutput->Branch("weight_ctwi_p5",&weight_ctwi_p5,"weight_ctwi_p5/F");

   //tOutput->Branch("weight_cbwi_m5",&weight_cbwi_m5,"weight_cbwi_m5/F");
   tOutput->Branch("weight_cbwi_m2",&weight_cbwi_m2,"weight_cbwi_m2/F");
   tOutput->Branch("weight_cbwi_m1",&weight_cbwi_m1,"weight_cbwi_m1/F");
   tOutput->Branch("weight_cbwi_p1",&weight_cbwi_p1,"weight_cbwi_p1/F");
   tOutput->Branch("weight_cbwi_p2",&weight_cbwi_p2,"weight_cbwi_p2/F");
   //tOutput->Branch("weight_cbwi_p5",&weight_cbwi_p5,"weight_cbwi_p5/F");

   


   if (fChain == 0) return;

	//Long64_t nentries = 100;
   Long64_t nentries = fChain->GetEntriesFast();

   TLorentzVector Ptop;
   TLorentzVector Pb;
   TLorentzVector Pw;
   TLorentzVector Pl;
   TLorentzVector Pnu;
   TLorentzVector Pqspec;
   double Pl_ID;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

	  //cout << "nParticles="<<Particle_<<endl;
	  for (int i=0; i<Particle_; i++){

      if (Rwgt_>0){
         cout<<"Rwgt_>0"<<endl;
         weight_SM = Rwgt_Weight[0];
         weight_ctwi_p1 = Rwgt_Weight[1];
         weight_ctwi_m1 = Rwgt_Weight[2];
         weight_ctwi_p2 = Rwgt_Weight[3];
         weight_ctwi_m2 = Rwgt_Weight[4];
         weight_cbwi_p1 = Rwgt_Weight[5];
         weight_cbwi_m1 = Rwgt_Weight[6];
         weight_cbwi_p2 = Rwgt_Weight[7];
         weight_cbwi_m2 = Rwgt_Weight[8];
      }


		  if (TMath::Abs(Particle_PID[i])==24 && Particle_Status[i]==2)
			   Pw.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);

		  if (TMath::Abs(Particle_PID[i])==5 && Particle_Status[i]==1)
			  Pb.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);

		  if (TMath::Abs(Particle_PID[i])<=5 && Particle_Status[i]==1)
			  Pqspec.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);

		  if ((TMath::Abs(Particle_PID[i])==11 || TMath::Abs(Particle_PID[i])==13) && Particle_Status[i]==1)
      {
        Pl.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);
        Pl_ID = abs(Particle_PID[i]);
        if (Pl_ID == 11)
        {
          nature_lepton = 1;
        }
        if (Pl_ID == 13)
        {
          nature_lepton = 2;
        }
      }


      if ((TMath::Abs(Particle_PID[i])==12 || TMath::Abs(Particle_PID[i])==14) && Particle_Status[i]==1)
  		  Pnu.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);
	  }

    Ptop = Pb + Pl + Pnu;

	  weight = Event_Weight[jentry];

	  /* SELECTION */

  	if ((Pl.Pt()<35 || TMath::Abs(Pl.Eta())>1.479) && Pl_ID==11) continue;

    if ((Pl.Pt()<26 || TMath::Abs(Pl.Eta())>2.4) && Pl_ID==13) continue;

	  if (Pqspec.Pt()<40 || TMath::Abs(Pqspec.Eta())>4.7) continue;
	  if (Pb.Pt()<40 || TMath::Abs(Pb.Eta())>2.4) continue;
    if (abs(Pqspec.Eta()) < 3.0 && 2.7 < abs(Pqspec.Eta()) && Pqspec.Pt() < 60) continue;


	  /* ANGLE RECONSTRUCTION */

	  TVector3 InvTopBoost;  InvTopBoost.SetXYZ(-Ptop.Px()/Ptop.E(),-Ptop.Py()/Ptop.E(),-Ptop.Pz()/Ptop.E());
	  TVector3 InvWBoost; InvWBoost.SetXYZ(-Pw.Px()/Pw.E(),-Pw.Py()/Pw.E(),-Pw.Pz()/Pw.E());

	  Pw.Boost(InvTopBoost);
	  Pb.Boost(InvTopBoost);
	  Pqspec.Boost(InvTopBoost);
	  //cout << "Top rest frame, W: Px="<<Pw.Px()<<" Py="<<Pw.Py()<<" Pz="<<Pw.Pz()<<endl;
	  //cout << "Top rest frame, b: Px="<<Pb.Px()<<" Py="<<Pb.Py()<<" Pz="<<Pb.Pz()<<endl;
	  //cout << "Top rest frame, q spec: Px="<<Pqspec.Px()<<" Py="<<Pqspec.Py()<<" Pz="<<Pqspec.Pz()<<endl;

	  TVector3 Zdir = Pw.Vect().Unit();
	  TVector3 PqspecUnit = Pqspec.Vect().Unit();
	  TVector3 Ydir = PqspecUnit.Cross(Zdir).Unit();
	  TVector3 Xdir = Ydir.Cross(Zdir);

	  sinTheta = PqspecUnit.Cross(Zdir).Mag();
	  cosTheta = PqspecUnit.Dot(Zdir);
	  if (sinTheta<0) sinTheta=-sinTheta;
	  //cout << "Angle Theta: cosTheta="<<cosTheta<<" sinTheta="<<sinTheta<<endl;

	  //cout << "Top rest frame, Xdir: X="<<Xdir.X()<<" Y="<<Xdir.Y()<<" Z="<<Xdir.Z()<< " Mag="<< Xdir.Mag() <<endl;
	  //cout << "Top rest frame, Ydir: X="<<Ydir.X()<<" Y="<<Ydir.Y()<<" Z="<<Ydir.Z()<< " Mag="<< Ydir.Mag() <<endl;
	  //cout << "Top rest frame, Zdir: X="<<Zdir.X()<<" Y="<<Zdir.Y()<<" Z="<<Zdir.Z()<< " Mag="<< Zdir.Mag() <<endl;

	  Pl.Boost(InvWBoost);
	  lepton_E_Wframe = Pl.E();

	  TVector3 PlUnit = Pl.Vect().Unit();
	  cosThetaStar = PlUnit.Dot(Zdir);
	  sinThetaStar = PlUnit.Cross(Zdir).Mag();
	  if (sinThetaStar<0) sinThetaStar=-sinThetaStar;
	  //cout << "Angle ThetaStar: cosThetaStar="<<cosThetaStar<<" sinThetaStar="<<sinThetaStar<<endl;

	  TVector3 PlUnit_PlaneXY = (PlUnit - (PlUnit.Dot(Zdir))*Zdir).Unit();
	  cosPhiStar = PlUnit_PlaneXY.Dot(Xdir);
	  sinPhiStar = PlUnit_PlaneXY.Dot(Ydir);
      //PhiStar;// = TMath::ACos(cosPhiStar);
	  if (sinPhiStar>0) PhiStar = TMath::ACos(cosPhiStar);
	  if (sinPhiStar<0) PhiStar = 2*TMath::Pi()-TMath::ACos(cosPhiStar);
	  //cout << "Angle PhiStar: cosPhiStar="<<cosPhiStar<<" sinPhiStar="<<sinPhiStar<<" PhiStar="<<PhiStar<<endl;

    top_pt = Ptop.Pt();
    W_pt = Pw.Pt();
    lepton_pt = Pl.Pt();
    top_mass = Ptop.M();
    W_mass = Pw.M();
    W_transverse_mass = Pw.Mt();

	  tOutput->Fill();
   }

   tOutput->Write();
   fOutput->Close();
}
