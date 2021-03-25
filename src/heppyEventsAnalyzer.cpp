#define heppyEventsAnalyzer_cxx
#include "../include/heppyEventsAnalyzer.hpp"
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

void heppyEventsAnalyzer::Loop()
{
  double w_mass = 80.4;
  float final_jet1_pt;
  float final_jet2_pt;
  float final_jet3_pt;
  float final_jet1_eta;
  float final_jet2_eta;
  float final_jet3_eta;

  float final_b1_pt;
  float final_b2_pt;
  float final_b3_pt;
  float final_b1_eta;
  float final_b2_eta;
  float final_b3_eta;

  float final_elec_iso;
  float final_muon_iso;
  float final_elec_pt;
  float final_muon_pt;
  float final_elec_eta;
  float final_muon_eta;

  float final_natureLepton;
  float final_region;
  float weight_one_lepton;
  float charge_lepton;

  double nb_elec_region1 = 0;
  double nb_muon_region1 = 0;
  double nb_antielec_region1 = 0;
  double nb_antimuon_region1 = 0;

  float M_T_top;
  float M_T_W;
  float M_T_tot;
  float jet_not_b_eta;

  float Delta_R;

  TLorentzVector Pl;
  float pt_lepton;
  float phi_lepton;
  float eta_lepton;
  float energy_lepton;
  float pz_lepton;

  TLorentzVector Pw;
  float reco_W_energy;
  float reco_W_pt;
  float reco_W_phi;
  float reco_W_eta;
  float reco_W_pz;

  TLorentzVector Pnu;
  float reco_nu_pz;
  float reco_nu_pt;
  float reco_nu_energy;
  float reco_nu_phi;
  float reco_nu_eta;

  TLorentzVector Ptop;
  float reco_top_pt;
  float reco_top_energy;
  float reco_top_phi;
  float reco_top_eta;
  float reco_top_pz;

  TLorentzVector Pb;
  TLorentzVector Pqspec;

  float boosted_reco_sinTheta, boosted_reco_cosTheta;
  float boosted_reco_sinThetaStar, boosted_reco_cosThetaStar;
  float boosted_reco_sinPhiStar, boosted_reco_cosPhiStar, boosted_reco_PhiStar;
  float boosted_reco_lepton_E_Wframe;
  float boosted_reco_W_pt, boosted_reco_lepton_pt;
  float boosted_reco_top_mass, boosted_reco_W_mass, boosted_reco_W_transverse_mass;

  int erreur = 0;
  TFile* fOutput = new TFile("output.root","RECREATE");
  TTree* tOutput = new TTree("events","events");

  tOutput->Branch("final_jet1_pt", &final_jet1_pt, "final_jet1_pt/F");
  tOutput->Branch("final_jet2_pt", &final_jet2_pt, "final_jet2_pt/F");
  tOutput->Branch("final_jet3_pt", &final_jet3_pt, "final_jet3_pt/F");

  tOutput->Branch("final_jet1_eta", &final_jet1_eta, "final_jet1_eta/F");
  tOutput->Branch("final_jet2_eta", &final_jet2_eta, "final_jet2_eta/F");
  tOutput->Branch("final_jet3_eta", &final_jet3_eta, "final_jet3_eta/F");

  tOutput->Branch("final_b1_pt", &final_b1_pt, "final_b1_pt/F");
  tOutput->Branch("final_b2_pt", &final_b2_pt, "final_b2_pt/F");
  tOutput->Branch("final_b3_pt", &final_b3_pt, "final_b3_pt/F");

  tOutput->Branch("final_b1_eta", &final_b1_eta, "final_b1_eta/F");
  tOutput->Branch("final_b2_eta", &final_b2_eta, "final_b2_eta/F");
  tOutput->Branch("final_b3_eta", &final_b3_eta, "final_b3_eta/F");

  tOutput->Branch("final_elec_iso", &final_elec_iso, "final_elec_iso/F");
  tOutput->Branch("final_muon_iso", &final_muon_iso, "final_muon_iso/F");
  tOutput->Branch("final_elec_pt", &final_elec_pt, "final_elec_pt/F");
  tOutput->Branch("final_muon_pt", &final_muon_pt, "final_muon_pt/F");
  tOutput->Branch("final_elec_eta", &final_elec_eta, "final_elec_eta/F");
  tOutput->Branch("final_muon_eta", &final_muon_eta, "final_muon_eta/F");

  tOutput->Branch("final_natureLepton", &final_natureLepton, "final_natureLepton/F");
  tOutput->Branch("final_region", &final_region, "final_region/F");
  tOutput->Branch("weight_one_lepton", &weight_one_lepton, "weight_one_lepton/F");
  tOutput->Branch("charge_lepton", &charge_lepton, "charge_lepton/F");

  tOutput->Branch("M_T_W", &M_T_W, "M_T_W/F");
  tOutput->Branch("M_T_top", &M_T_top, "M_T_top/F");
  tOutput->Branch("M_T_tot", &M_T_tot, "M_T_tot/F");
  tOutput->Branch("jet_not_b_eta", &jet_not_b_eta, "jet_not_b_eta/F");
  tOutput->Branch("Delta_R", &Delta_R, "Delta_R/F");

  tOutput->Branch("pt_lepton", &pt_lepton, "pt_lepton/F");
  tOutput->Branch("phi_lepton", &phi_lepton, "phi_lepton/F");
  tOutput->Branch("eta_lepton", &eta_lepton, "eta_lepton/F");
  tOutput->Branch("energy_lepton", &energy_lepton, "energy_lepton/F");
  tOutput->Branch("pz_lepton", &pz_lepton, "pz_lepton/F");



  tOutput->Branch("reco_nu_energy", &reco_nu_energy, "reco_nu_energy/F");
  tOutput->Branch("reco_nu_pt", &reco_nu_pt, "reco_nu_pt/F");
  tOutput->Branch("reco_nu_phi", &reco_nu_phi, "reco_nu_phi/F");
  tOutput->Branch("reco_nu_eta", &reco_nu_eta, "reco_nu_eta/F");
  tOutput->Branch("reco_nu_pz", &reco_nu_pz, "reco_nu_pz/F");

  tOutput->Branch("reco_W_energy", &reco_W_energy, "reco_W_energy/F");
  tOutput->Branch("reco_W_pt", &reco_W_pt, "reco_W_pt/F");
  tOutput->Branch("reco_W_phi", &reco_W_phi, "reco_W_phi/F");
  tOutput->Branch("reco_W_eta", &reco_W_eta, "reco_W_eta/F");
  tOutput->Branch("reco_W_pz", &reco_W_pz, "reco_W_pz/F");

  tOutput->Branch("reco_top_energy", &reco_top_energy, "reco_top_energy/F");
  tOutput->Branch("reco_top_pt", &reco_top_pt, "reco_top_pt/F");
  tOutput->Branch("reco_top_phi", &reco_top_phi, "reco_top_phi/F");
  tOutput->Branch("reco_top_eta", &reco_top_eta, "reco_top_eta/F");
  tOutput->Branch("reco_top_pz", &reco_top_pz, "reco_top_pz/F");



  tOutput->Branch("boosted_reco_W_pt", &boosted_reco_W_pt, "boosted_reco_W_pt/F");
  tOutput->Branch("boosted_reco_lepton_pt", &boosted_reco_lepton_pt, "boosted_reco_lepton_pt/F");
  tOutput->Branch("boosted_reco_top_mass", &boosted_reco_top_mass, "boosted_reco_top_mass/F");
  tOutput->Branch("boosted_reco_W_mass", &boosted_reco_W_mass, "boosted_reco_W_mass/F");
  tOutput->Branch("boosted_reco_W_transverse_mass", &boosted_reco_W_transverse_mass, "boosted_reco_W_transverse_mass/F");
  tOutput->Branch("boosted_reco_lepton_E_Wframe", &boosted_reco_lepton_E_Wframe, "boosted_reco_lepton_E_Wframe/F");
  tOutput->Branch("boosted_reco_PhiStar", &boosted_reco_PhiStar, "boosted_reco_PhiStar/F");
  tOutput->Branch("boosted_reco_cosPhiStar", &boosted_reco_cosPhiStar, "boosted_reco_cosPhiStar/F");
  tOutput->Branch("boosted_reco_sinPhiStar", &boosted_reco_sinPhiStar, "boosted_reco_sinPhiStar/F");
  tOutput->Branch("boosted_reco_cosThetaStar", &boosted_reco_cosThetaStar, "boosted_reco_cosThetaStar/F");
  tOutput->Branch("boosted_reco_sinThetaStar", &boosted_reco_sinThetaStar, "boosted_reco_sinThetaStar/F");
  tOutput->Branch("boosted_reco_cosTheta", &boosted_reco_cosTheta, "boosted_reco_cosTheta/F");
  tOutput->Branch("boosted_reco_sinTheta", &boosted_reco_sinTheta, "boosted_reco_sinTheta/F");





   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      final_jet1_pt = -10000;
      final_jet2_pt = -10000;
      final_jet3_pt = -10000;
      final_jet1_eta = -10000;
      final_jet2_eta = -10000;
      final_jet3_eta = -10000;

      final_b1_pt = -10000;
      final_b2_pt = -10000;
      final_b3_pt = -10000;
      final_b1_eta = -10000;
      final_b2_eta = -10000;
      final_b3_eta = -10000;

      final_elec_iso = -10000;
      final_muon_iso = -10000;
      final_elec_pt = -10000;
      final_muon_pt = -10000;
      final_elec_eta = -10000;
      final_muon_eta = -10000;

      final_natureLepton = -10000;
      final_region = -10000;
      weight_one_lepton = -10000;
      charge_lepton = -10000;
      final_natureLepton = natureLepton;

      pt_lepton = -10000;
      phi_lepton = -10000;
      eta_lepton = -10000;
      energy_lepton = -10000;

      pt_lepton = -10000;
      phi_lepton = -10000;
      eta_lepton = -10000;
      energy_lepton = -10000;

      reco_W_energy = -10000;
      reco_W_pt = -10000;
      reco_W_phi = -10000;
      reco_W_eta = -10000;
      reco_W_pz = -10000;

      reco_nu_pz = -10000;
      reco_nu_pt = -10000;
      reco_nu_energy = -10000;
      reco_nu_phi = -10000;
      reco_nu_eta = -10000;

      reco_top_pt = -10000;
      reco_top_energy = -10000;
      reco_top_phi = -10000;
      reco_top_eta = -10000;
      reco_top_pz = -10000;

      boosted_reco_sinTheta = -10000;
      boosted_reco_cosTheta = -10000;
      boosted_reco_sinThetaStar  = -10000;
      boosted_reco_cosThetaStar = -10000;
      boosted_reco_sinPhiStar  = -10000;
      boosted_reco_cosPhiStar = -10000;
      boosted_reco_PhiStar = -10000;
      boosted_reco_lepton_E_Wframe  = -10000;
      boosted_reco_W_pt = -10000;
      boosted_reco_lepton_pt = -10000;
      boosted_reco_top_mass  = -10000;
      boosted_reco_W_mass = -10000;
      boosted_reco_W_transverse_mass = -10000;

        if ((2.7 < abs(j1_eta) && abs(j1_eta) < 3.0 && j1_pt < 60 ) || (2.7 < abs(j2_eta) && abs(j2_eta)< 3.0 && j2_pt < 60 ) || (2.7 < abs(j3_eta) && abs(j3_eta) <3.0 && j3_pt < 60 )  )
        {
          continue;
        }

        if (is_data == 0)
        {
          weight_one_lepton = weight*weight_generator;
        }
        else {weight_one_lepton = 1;}


        final_region = region;
        final_jet1_pt = j1_pt;
        final_jet2_pt = j2_pt;
        final_jet3_pt = j3_pt;
        final_jet1_eta = j1_eta;
        final_jet2_eta = j2_eta;
        final_jet3_eta = j3_eta;

        final_b1_pt = b1_pt;
        final_b2_pt = b2_pt;
        final_b3_pt = b3_pt;
        final_b1_eta = b1_eta;
        final_b2_eta = b2_eta;
        final_b3_eta = b3_eta;

        final_elec_iso = iso_elec;
        final_muon_iso = iso_muon;
        final_elec_pt = pt_elec;
        final_muon_pt = pt_muon;
        final_elec_eta = eta_elec;
        final_muon_eta = eta_muon;

        if (j1_pt == b1_pt)
        {
          jet_not_b_eta = j2_eta;
        }
        if (j2_pt == b1_pt)
        {
          jet_not_b_eta = j1_eta;
        }

        if (region == 1)
          {

            if (natureLepton == 1)
            {

              if (q_elec == -1)
              {
                nb_elec_region1 += weight_one_lepton;

              }
              if (q_elec == 1)
              {
                nb_antielec_region1 += weight_one_lepton;
              }
            }
            if (natureLepton == 2)
            {

              if (q_muon == -1)
              {
                nb_muon_region1 += weight_one_lepton;
              }
              if (q_muon == 1)
              {
                nb_antimuon_region1 += weight_one_lepton;
              }
            }

          }

      //--------------------------W reconstruction-------------------------//

      if (j1_pt == b1_pt)
      {
        Pqspec.SetPxPyPzE(j2_pt*cos(j2_phi),j2_pt*sin(j2_phi), j2_pt*sinh(j2_eta), j2_energy);
      }
      if (j2_pt == b1_pt)
      {
        Pqspec.SetPxPyPzE(j1_pt*cos(j1_phi),j1_pt*sin(j1_phi), j1_pt*sinh(j1_eta), j1_energy);
      }

      if (natureLepton == 1)
      {
        Pl.SetPxPyPzE(pt_elec*cos(phi_elec),pt_elec*sin(phi_elec), pt_elec*sinh(eta_elec), energy_elec);
        charge_lepton = q_elec;

      }
      if (natureLepton == 2)
      {
        Pl.SetPxPyPzE(pt_muon*cos(phi_muon),pt_muon*sin(phi_muon), pt_muon*sinh(eta_muon), energy_muon);
        charge_lepton = q_muon;
      }

      M_T_W = sqrt(pow(met*cos(metphi)+Pl.Px(),2) + pow(met*sin(metphi)+Pl.Py(),2) );

      pt_lepton = Pl.Pt();
      phi_lepton = Pl.Phi();
      eta_lepton = Pl.Eta();
      energy_lepton = Pl.Energy();
      pz_lepton = Pl.Pz();

      double lambda = w_mass*w_mass/2 + met*cos(metphi)*Pl.Px()+met*sin(metphi)*Pl.Py();

      double Delta = 4*pow(lambda*Pl.Pz(),2) - 4*Pl.Pt()*Pl.Pt()*(Pl.Energy()*Pl.Energy()*met*met-lambda*lambda);


    //-------------Delta > 0------------//
  if (Delta > 0)
    {
        //---------Compute reco_nu_pt---------//

        reco_nu_pt = met;
        reco_nu_phi = metphi;

        double reco_nu_pz_up = lambda*Pl.Pz()/(Pl.Pt()*Pl.Pt()) + sqrt((lambda*lambda*Pl.Pz()*Pl.Pz())/(Pl.Pt()*Pl.Pt()*Pl.Pt()*Pl.Pt()) - (Pl.Energy()*Pl.Energy()*met*met-lambda*lambda)/(Pl.Pt()*Pl.Pt()));
        double reco_nu_pz_down = lambda*Pl.Pz()/(Pl.Pt()*Pl.Pt()) - sqrt(pow(lambda*Pl.Pz()/(Pl.Pt()*Pl.Pt()),2) - (Pl.Energy()*Pl.Energy()*met*met-lambda*lambda)/(Pl.Pt()*Pl.Pt()));
        //------------Compute reco_nu_pz--------------//

          if (abs(reco_nu_pz_up) > abs(reco_nu_pz_down))
            {
              reco_nu_pz = reco_nu_pz_down;
            }
          else
            {
              reco_nu_pz = reco_nu_pz_up;
            }

          reco_nu_energy = sqrt(met*met+reco_nu_pz*reco_nu_pz);

          Pnu.SetPxPyPzE(reco_nu_pt*cos(reco_nu_phi),reco_nu_pt*sin(reco_nu_phi), reco_nu_pz, reco_nu_energy);

      }
    //----------------Delta <= 0 -----------//
   if (Delta < 0)
  {
      //------------Compute reco_nu_pt-----------//

      double reco_nu_pt_up = sqrt(2)*abs(w_mass + Pl.Pt()/sqrt(2));
      double reco_nu_pt_down = sqrt(2)*abs(w_mass - Pl.Pt()/sqrt(2));

      if (reco_nu_pt_down < 0)
      {
        reco_nu_pt = reco_nu_pt_up;
      }
      else if (abs(met - reco_nu_pt_down) < abs(met - reco_nu_pt_up))
      {
        reco_nu_pt = reco_nu_pt_down;
      }
      else if (abs(met - reco_nu_pt_down) > abs(met - reco_nu_pt_up))
      {
        reco_nu_pt = reco_nu_pt_up;
      }


      //-------Compute phi of nu---------//

      double cos_theta_lep_nu = (w_mass*w_mass - reco_nu_pt*reco_nu_pt - (Pl.Pt()*Pl.Pt()))/(2*Pl.Pt()*reco_nu_pt);



      if (abs(cos_theta_lep_nu)>1 && abs(met - reco_nu_pt_down) < abs(met - reco_nu_pt_up) && reco_nu_pt_up > 0 )
      {
        reco_nu_pt = reco_nu_pt_up;
      }
      if (abs(cos_theta_lep_nu)>1 && abs(met - reco_nu_pt_down) > abs(met - reco_nu_pt_up) && reco_nu_pt_down > 0)
      {
        reco_nu_pt = reco_nu_pt_down;
      }

      cos_theta_lep_nu = (w_mass*w_mass - reco_nu_pt*reco_nu_pt - Pl.Pt()*Pl.Pt())/(2*Pl.Pt()*reco_nu_pt);

      if (cos_theta_lep_nu >=0) reco_nu_phi = phi_lepton + acos(cos_theta_lep_nu);

      if (cos_theta_lep_nu<0) reco_nu_phi = phi_lepton + (2*M_PI - acos(cos_theta_lep_nu));

      if (isnan(reco_nu_phi)) continue;


      //--------Compute reco_nu_pz when met != reco_nu_pt-------//
      lambda = w_mass*w_mass/2 + reco_nu_pt*Pl.Pt()*cos_theta_lep_nu;

        reco_nu_pz = lambda*Pl.Pz()/(Pl.Pt()*Pl.Pt());
        reco_nu_energy = sqrt(reco_nu_pt*reco_nu_pt + reco_nu_pz*reco_nu_pz);


        Pnu.SetPxPyPzE(reco_nu_pt*cos(reco_nu_phi),reco_nu_pt*sin(reco_nu_phi), reco_nu_pz, reco_nu_energy);

   }

    reco_nu_pt = Pnu.Pt();
    reco_nu_phi = Pnu.Phi();
    reco_nu_eta = Pnu.Eta();
    reco_nu_pz = Pnu.Pz();
    reco_nu_energy = Pnu.Energy();


    Pw = Pl + Pnu;
    reco_W_pt = Pw.Pt();
    reco_W_phi = Pw.Phi();
    reco_W_eta = Pw.Eta();
    reco_W_pz = Pw.Pz();
    reco_W_energy = Pw.Energy();



    //--------------Top reconstruction-----------//
    Pb.SetPtEtaPhiE(b1_pt, b1_eta, b1_phi, b1_energy);

    Ptop = Pw + Pb;
    reco_top_pt = Ptop.Pt();
    reco_top_pz = Ptop.Pz();
    reco_top_eta = Ptop.Eta();
    reco_top_phi = Ptop.Phi();
    reco_top_energy = Ptop.Energy();

    if((abs(Ptop.Pt())> 100000 || abs(Ptop.Phi())> 100000  || abs(Ptop.Eta())> 100000 || abs(Ptop.Pz())> 100000 || abs(Ptop.Energy())> 100000) && final_region == 1)
    {
      erreur++;
    }

    M_T_top = sqrt(pow(met*cos(metphi)+Pl.Px()+Pb.Px(),2) + pow(met*sin(metphi)+Pl.Py()+Pb.Py(),2) );
    M_T_tot = sqrt(pow(met*cos(metphi)+Pl.Px()+j1_pt*cos(j1_phi)+j2_pt*cos(j2_phi)+j3_pt*cos(j3_phi),2) + pow(met*sin(metphi)+Pl.Py()+j1_pt*sin(j1_phi)+j2_pt*sin(j2_phi)+j3_pt*sin(j3_phi),2) );

    if (region == 1)
    {
      Delta_R = Pl.DeltaR(Pqspec);
    }


    /*if (j1_pt == b1_pt && j2_pt == b2_pt) most_forward_not_b_eta = j3_eta;

    if (j2_pt == b1_pt && j3_pt == b2_pt) most_forward_not_b_eta = j1_eta;*/










    TVector3 InvTopBoost;  InvTopBoost.SetXYZ(-Ptop.Px()/Ptop.E(),-Ptop.Py()/Ptop.E(),-Ptop.Pz()/Ptop.E());
    TVector3 InvWBoost; InvWBoost.SetXYZ(-Pw.Px()/Pw.E(),-Pw.Py()/Pw.E(),-Pw.Pz()/Pw.E());

    Pw.Boost(InvTopBoost);
    Pb.Boost(InvTopBoost);
    Pqspec.Boost(InvTopBoost);

    TVector3 Zdir = Pw.Vect().Unit();
    TVector3 PqspecUnit = Pqspec.Vect().Unit();
    TVector3 Ydir = PqspecUnit.Cross(Zdir).Unit();
    TVector3 Xdir = Ydir.Cross(Zdir);

    boosted_reco_sinTheta = PqspecUnit.Cross(Zdir).Mag();
    boosted_reco_cosTheta = PqspecUnit.Dot(Zdir);
    if (boosted_reco_sinTheta<0) boosted_reco_sinTheta=-boosted_reco_sinTheta;

    Pl.Boost(InvWBoost);
    boosted_reco_lepton_E_Wframe = Pl.E();

    TVector3 PlUnit = Pl.Vect().Unit();
    boosted_reco_cosThetaStar = PlUnit.Dot(Zdir);
    boosted_reco_sinThetaStar = PlUnit.Cross(Zdir).Mag();

    if (boosted_reco_sinThetaStar<0) boosted_reco_sinThetaStar=-boosted_reco_sinThetaStar;

    TVector3 PlUnit_PlaneXY = (PlUnit - (PlUnit.Dot(Zdir))*Zdir).Unit();
    boosted_reco_cosPhiStar = PlUnit_PlaneXY.Dot(Xdir);
    boosted_reco_sinPhiStar = PlUnit_PlaneXY.Dot(Ydir);

    //boosted_reco_PhiStar = TMath::ACos(boosted_reco_cosPhiStar);
    //boosted_reco_PhiStar = 2*TMath::Pi()-TMath::ACos(boosted_reco_cosPhiStar);
    if (boosted_reco_sinPhiStar>=0) boosted_reco_PhiStar = TMath::ACos(boosted_reco_cosPhiStar);
    if (boosted_reco_sinPhiStar<0) boosted_reco_PhiStar = 2*TMath::Pi()-TMath::ACos(boosted_reco_cosPhiStar);
    //cout<<"sinPhiStar ="<<boosted_reco_sinPhiStar<<endl;

    boosted_reco_W_pt = Pw.Pt();
    boosted_reco_lepton_pt = Pl.Pt();
    boosted_reco_top_mass = Ptop.M();
    boosted_reco_W_mass = Pw.M();
    boosted_reco_W_transverse_mass = Pw.Mt();


  tOutput->Fill();
   }
   //cout<<"erreur = "<<erreur<<endl;
      /*cout<<"Number of electrons in W + jets control region :"<<nb_elec_region0<<endl;
   cout<<"Number of muons in W + jets control region :"<<nb_muon_region0<<endl;
   cout<<"Number of positrons in W + jets control region :"<<nb_antielec_region0<<endl;
   cout<<"Number of antimuons in W + jets control region :"<<nb_antimuon_region0<<endl;*/
   /*cout<<"Number of electron in single top signal region :"<<nb_elec_region1<<endl;
   cout<<"Number of positrons in single top signal region :"<<nb_antielec_region1<<endl;
   cout<<"Number of muons in single top signal region :"<<nb_muon_region1<<endl;
   cout<<"Number of antimuons in single top signal region :"<<nb_antimuon_region1<<endl;*/
   /*cout<<"Number of electron in ttbar control region :"<<nb_elec_region2<<endl;
   cout<<"Number of muons in ttbar control region :"<<nb_muon_region2<<endl;
   cout<<"Number of positrons in ttbar control region :"<<nb_antielec_region2<<endl;
   cout<<"Number of antimuons in ttbar control region :"<<nb_antimuon_region2<<endl;*/
   tOutput->Write();
   fOutput->Close();

}
