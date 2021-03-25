#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <sstream>
#include "TStyle.h"
#include <iostream>
#include <TStyle.h>
#include <THStack.h>
#include <TFractionFitter.h>
#include <TF1.h>

using namespace std;

TH1D* GetHistoWeight(TTree* t, string variable, int nbins, double xmin, double xmax, string cut, string name)
{
        string sxmin, sxmax, snbins;
        stringstream ss[3];

        ss[0] << xmin;
        ss[0] >> sxmin;
        ss[1] << xmax;
        ss[1] >> sxmax;
        ss[2] << nbins;
        ss[2] >> snbins;
        string variablenew = variable + " >> h(" + snbins + "," + sxmin + "," + sxmax + ")";


        string cutnew = "weight_one_lepton*(" + cut + ")";
        t->Draw(variablenew.c_str(), cutnew.c_str());
        TH1D *histo = (TH1D*)gDirectory->Get("h");
        histo->SetBinContent(0,0);
		if (histo->GetEntries()==0) return histo;

		double underflow = histo->GetBinContent(0);
		//cout << "underflow="<<underflow<<endl;
		double val = 0;
		if (underflow>0) {
			val = histo->GetBinContent(1);
			//histo->SetBinContent(1, val+underflow);
		}
		double overflow = histo->GetBinContent(nbins+1);
		if (overflow>0) {
		  val = histo->GetBinContent(nbins);
		  histo->SetBinContent(nbins+1, 0);
		  histo->SetBinContent(nbins, val+overflow);

		}

		//cout << "Area="<<histo->Integral()<<endl;
		//cout << "Nevents="<<histo->GetEntries()<<endl;
        histo->SetName(name.c_str());
        histo->SetTitle(name.c_str());
        return histo;
}

double normalizeWithChi(TTree* t1, TTree* t2, TTree* t3, TTree* t4 , TTree* t5, TTree* t6, TTree* t7, TTree* t8, TTree* t9, TTree* t10, TTree* t11, TTree* t12,
   string variable, int nbins, double xmin, double xmax, string selection, double param0_elec, double param1_elec, double param0_muon, double param1_muon, double syst_lumi, double syst_singletop, double syst_ttbar_tw, double syst_w_z_gamma, double syst_QCD, string signal, string legendX, string legendY, string Name, string region)
{

  double lumi = 41500;
  TH1D* Histo_1_elec = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 1", "");
  TH1D* Histo_2_elec = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 1", "");
  TH1D* Histo_3_elec = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 1", "");
  TH1D* Histo_4_elec = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 1", "");
  TH1D* Histo_5_elec = GetHistoWeight(t5, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 1", "");
  TH1D* Histo_6_elec = GetHistoWeight(t6, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 1", "");
  TH1D* Histo_7_elec = GetHistoWeight(t7, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 1", "");
  TH1D* Histo_8_elec = GetHistoWeight(t8, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 1", "");
  TH1D* Histo_9_elec = GetHistoWeight(t9, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 1", "");
  TH1D* Histo_10_elec = GetHistoWeight(t10, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 1", "");
  TH1D* DATA_elec = GetHistoWeight(t12, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 1", "");

  TH1D* Histo_1_muon = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 2", "");
  TH1D* Histo_2_muon = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 2", "");
  TH1D* Histo_3_muon = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 2", "");
  TH1D* Histo_4_muon = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 2", "");
  TH1D* Histo_5_muon = GetHistoWeight(t5, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 2", "");
  TH1D* Histo_6_muon = GetHistoWeight(t6, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 2", "");
  TH1D* Histo_7_muon = GetHistoWeight(t7, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 2", "");
  TH1D* Histo_8_muon = GetHistoWeight(t8, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 2", "");
  TH1D* Histo_9_muon = GetHistoWeight(t9, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 2", "");
  TH1D* Histo_10_muon = GetHistoWeight(t11, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 2", "");
  TH1D* DATA_muon = GetHistoWeight(t12, variable, nbins, xmin, xmax, selection+"&& final_natureLepton == 2", "");

  double events_MC_t_top = 5982064;
  double events_MC_t_antitop = 3675910;
  double events_MC_DY_50 = 491250940213;
  double events_MC_DY_1050 = 315866;
  double events_MC_semilep = 32414370098;
  double events_MC_tW_antitop = 270762750.172;
  double events_MC_tW_top = 277241050.839;
  double events_MC_W_jets = 29981320;
  double events_MC_dilep = 4984888995.13;

  double cross_MC_t_top = 136.02;
  double cross_MC_t_antitop = 80.95;
  double cross_MC_DY_50 = 6025.2;
  double cross_MC_DY_1050 = 22635.1;
  double cross_MC_semilep = 365.3;
  double cross_MC_tW_antitop = 35.9;
  double cross_MC_tW_top = 35.9;
  double cross_MC_W_jets = 61526.7;
  double cross_MC_dilep = 88.2;

  Histo_1_elec->Scale(1.0/events_MC_t_top*cross_MC_t_top*lumi);
  Histo_2_elec->Scale(1.0/events_MC_t_antitop*cross_MC_t_antitop*lumi);
  Histo_3_elec->Scale(1.0/events_MC_DY_50*cross_MC_DY_50*lumi);
  Histo_4_elec->Scale(1.0/events_MC_DY_1050*cross_MC_DY_1050*lumi);
  Histo_5_elec->Scale(1.0/events_MC_semilep*cross_MC_semilep*lumi);
  Histo_6_elec->Scale(1.0/events_MC_tW_antitop*cross_MC_tW_antitop*lumi);
  Histo_7_elec->Scale(1.0/events_MC_tW_top*cross_MC_tW_top*lumi);
  Histo_8_elec->Scale(1.0/events_MC_W_jets*cross_MC_W_jets*lumi);
  Histo_9_elec->Scale(1.0/events_MC_dilep*cross_MC_dilep*lumi);

  Histo_1_muon->Scale(1.0/events_MC_t_top*cross_MC_t_top*lumi);
  Histo_2_muon->Scale(1.0/events_MC_t_antitop*cross_MC_t_antitop*lumi);
  Histo_3_muon->Scale(1.0/events_MC_DY_50*cross_MC_DY_50*lumi);
  Histo_4_muon->Scale(1.0/events_MC_DY_1050*cross_MC_DY_1050*lumi);
  Histo_5_muon->Scale(1.0/events_MC_semilep*cross_MC_semilep*lumi);
  Histo_6_muon->Scale(1.0/events_MC_tW_antitop*cross_MC_tW_antitop*lumi);
  Histo_7_muon->Scale(1.0/events_MC_tW_top*cross_MC_tW_top*lumi);
  Histo_8_muon->Scale(1.0/events_MC_W_jets*cross_MC_W_jets*lumi);
  Histo_9_muon->Scale(1.0/events_MC_dilep*cross_MC_dilep*lumi);

  TH1D *singletop_elec = new TH1D ("stack single top", "", nbins, xmin, xmax );
  TH1D *ttbar_tw_elec= new TH1D ("stack ttbar_tw", "", nbins, xmin, xmax );
  TH1D *QCD_elec = new TH1D ("stack QCD", "", nbins, xmin, xmax );
  TH1D *w_z_gamma_elec = new TH1D ("stack w_z_gamma", "", nbins, xmin, xmax );
  TH1D *all_others_elec = new TH1D ("all_others", "", nbins, xmin, xmax );
  TH1D *ASIMOV_elec = new TH1D ("ASIMOV", "", nbins, xmin, xmax );
  TH1D *EFT_signal = new TH1D ("EFT_signal", "", nbins, xmin, xmax );

  TH1D *singletop_muon = new TH1D ("stack single top", "", nbins, xmin, xmax );
  TH1D *ttbar_tw_muon= new TH1D ("stack ttbar_tw", "", nbins, xmin, xmax );
  TH1D *QCD_muon = new TH1D ("stack QCD", "", nbins, xmin, xmax );
  TH1D *w_z_gamma_muon = new TH1D ("stack w_z_gamma", "", nbins, xmin, xmax );
  TH1D *all_others_muon = new TH1D ("all_others", "", nbins, xmin, xmax );
  TH1D *ASIMOV_muon = new TH1D ("ASIMOV", "", nbins, xmin, xmax );



  singletop_elec->Add(Histo_1_elec);
  singletop_elec->Add(Histo_2_elec);
  ttbar_tw_elec->Add(Histo_5_elec);
  ttbar_tw_elec->Add(Histo_9_elec);
  ttbar_tw_elec->Add(Histo_6_elec);
  ttbar_tw_elec->Add(Histo_7_elec);
  w_z_gamma_elec->Add(Histo_3_elec);
  w_z_gamma_elec->Add(Histo_4_elec);
  w_z_gamma_elec->Add(Histo_8_elec);
  QCD_elec->Add(Histo_10_elec);
  all_others_elec->Add(singletop_elec);
  all_others_elec->Add(ttbar_tw_elec);
  all_others_elec->Add(w_z_gamma_elec);

  singletop_muon->Add(Histo_1_muon);
  singletop_muon->Add(Histo_2_muon);
  ttbar_tw_muon->Add(Histo_5_muon);
  ttbar_tw_muon->Add(Histo_9_muon);
  ttbar_tw_muon->Add(Histo_6_muon);
  ttbar_tw_muon->Add(Histo_7_muon);
  w_z_gamma_muon->Add(Histo_3_muon);
  w_z_gamma_muon->Add(Histo_4_muon);
  w_z_gamma_muon->Add(Histo_8_muon);
  QCD_muon->Add(Histo_10_muon);
  all_others_muon->Add(singletop_muon);
  all_others_muon->Add(ttbar_tw_muon);
  all_others_muon->Add(w_z_gamma_muon);


  double integ_all_others_elec = all_others_elec->Integral(1,nbins);
  double integ_QCD_elec = QCD_elec->Integral(1,nbins);
  double integ_data_elec = DATA_elec->Integral(1,nbins);

  double integ_all_others_muon = all_others_muon->Integral(1,nbins);
  double integ_QCD_muon = QCD_muon->Integral(1,nbins);
  double integ_data_muon = DATA_muon->Integral(1,nbins);

  double newscale0_elec = integ_data_elec*param0_elec/integ_all_others_elec;
  double newscale1_elec = integ_data_elec*param1_elec/integ_QCD_elec;

  double newscale0_muon = integ_data_muon*param0_muon/integ_all_others_muon;
  double newscale1_muon = integ_data_muon*param1_muon/integ_QCD_muon;

  singletop_elec->Scale(newscale0_elec);
  ttbar_tw_elec->Scale(newscale0_elec);
  w_z_gamma_elec->Scale(newscale0_elec);
  QCD_elec->Scale(newscale1_elec);

  ASIMOV_elec->Add(singletop_elec);
  ASIMOV_elec->Add(ttbar_tw_elec);
  ASIMOV_elec->Add(w_z_gamma_elec);
  ASIMOV_elec->Add(QCD_elec);

  singletop_muon->Scale(newscale0_muon);
  ttbar_tw_muon->Scale(newscale0_muon);
  w_z_gamma_muon->Scale(newscale0_muon);
  QCD_muon->Scale(newscale1_muon);

  ASIMOV_muon->Add(singletop_muon);
  ASIMOV_muon->Add(ttbar_tw_muon);
  ASIMOV_muon->Add(w_z_gamma_muon);
  ASIMOV_muon->Add(QCD_muon);


  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  TFile* EFT_file_elec = TFile::Open(("./"+signal+"_elec_TF1.root").c_str());
  TH1D** EFT_histo_elec = new TH1D*[nbins];
  TF1** ratio_elec = new TF1*[nbins];


  TFile* EFT_file_muon = TFile::Open(("./"+signal+"_muon_TF1.root").c_str());
  TH1D** EFT_histo_muon = new TH1D*[nbins];
  TF1** ratio_muon = new TF1*[nbins];

  for(int i = 1 ; i <= nbins; i++)
  {
    ratio_elec[i-1] = (TF1*)EFT_file_elec->Get(("bin_content_par1_"+to_string(i)).c_str());

    EFT_histo_elec[i-1] = (TH1D*)ratio_elec[i-1]->CreateHistogram();

    ratio_muon[i-1] = (TF1*)EFT_file_muon->Get(("bin_content_par1_"+to_string(i)).c_str());
    EFT_histo_muon[i-1] = (TH1D*)ratio_muon[i-1]->CreateHistogram();
  }

  for (int p = 1; p <= nbins ; p++)
    {
      EFT_signal->SetBinContent(p, singletop_elec->GetBinContent(p)*ratio_elec[p-1]->Eval(2.0));
    }

  double error_Wjets_tot[nbins];
  double error_ttbar_tw_tot[nbins];
  double error_tchan_tot[nbins];
  double error_multijet_tot[nbins];

    for (int i = 0 ; i< nbins ; i++)
    {

      error_Wjets_tot[i] = sqrt(pow(3200,2) + pow(2800,2))/(33400+30700)*w_z_gamma_elec->GetBinContent(i+1);
      error_ttbar_tw_tot[i] = sqrt(pow(1400,2)+ pow(1500,2))/(84500+84800)*ttbar_tw_elec->GetBinContent(i+1);
      error_tchan_tot[i] = sqrt(pow(820,2)+pow(880,2)+pow(3,2)+pow(2,2))/(17720+27+25+11460)*singletop_elec->GetBinContent(i+1);
      error_multijet_tot[i] = 0.2*QCD_elec->GetBinContent(i+1);

    }

  TH1D *all_error_up = new TH1D("error_up","error_up",nbins,xmin,xmax);
  TH1D *all_histo = new TH1D ("multijet", "multijet", nbins, xmin, xmax );
  TH1D *all_m1_histo = new TH1D ("w_z_gamma", "w_z_gamma", nbins, xmin, xmax );
  TH1D *all_m2_histo = new TH1D ("ttbar_tw", "ttbar_tw", nbins, xmin, xmax );
  TH1D *all_m3_histo = new TH1D ("singletop", "singletop", nbins, xmin, xmax );

  all_histo->Add(singletop_elec);
  all_histo->Add(ttbar_tw_elec);
  all_histo->Add(w_z_gamma_elec);
  all_histo->Add(QCD_elec);

  all_m1_histo->Add(singletop_elec);
  all_m1_histo->Add(ttbar_tw_elec);
  all_m1_histo->Add(w_z_gamma_elec);

  all_m2_histo->Add(singletop_elec);
  all_m2_histo->Add(ttbar_tw_elec);

  all_m3_histo->Add(singletop_elec);

  all_histo->SetLineColor(kBlack);
  all_histo->SetFillColor(kGray);

  all_m1_histo->SetLineColor(kBlack);
  all_m1_histo->SetFillColor(8);

  all_m2_histo->SetLineColor(kBlack);
  all_m2_histo->SetFillColor(kOrange);

  all_m3_histo->SetLineColor(kBlack);
  all_m3_histo->SetFillColor(kRed);

  DATA_elec->SetLineWidth(2);
  DATA_elec->SetLineColor(kBlack);

  EFT_signal->SetLineWidth(3);
  EFT_signal->SetLineColor(kBlue);

  all_error_up->SetFillStyle(3005);
  all_error_up->SetFillColor(kBlack);
  all_error_up->SetMarkerStyle(0);

  for (int i = 0 ; i<nbins ; i++)
  {
    all_error_up->SetBinContent(i+1, all_histo->GetBinContent(i+1));
    all_error_up->SetBinError(i+1,sqrt( pow(error_Wjets_tot[i],2) + pow(error_ttbar_tw_tot[i],2) + pow(error_tchan_tot[i],2) +pow(error_multijet_tot[i],2)));
  }

  all_histo->SetXTitle(legendX.c_str());
  all_histo->SetYTitle(legendY.c_str());

  double max = all_histo->GetMaximum();
  all_histo->SetAxisRange(0,max*1.8,"Y");

  all_histo->Draw("HIST");
  all_error_up->Draw("E2 SAME");
  all_m1_histo->Draw("SAME HIST");
  all_m2_histo->Draw("SAME HIST");
  all_m3_histo->Draw("SAME HIST");
  EFT_signal->Draw("SAME");
  DATA_elec->Draw("E SAME");

  string legendtitle1 = region+" Region";
  float lx0 = 0.65;
  float ly0 = 0.65;
  float lx1 = 0.99;
  float ly1 = 0.99;

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle1.c_str());
  string legendEntry1 = "DATA";
  string legendEntry2 = "Multijet (DATA)";
  string legendEntry3 = "W/Z/#gamma* + jets";
  string legendEntry4 = "t#bar{t}/tw";
  string legendEntry5 = "Single top t-channel";
  string legendEntry6 = "C_{tW}^{I} = 2";



  legend->SetFillColor(kWhite);
  legend->AddEntry(DATA_elec->GetName(), legendEntry1.c_str(), "l");
  legend->AddEntry(all_histo->GetName(), legendEntry2.c_str(), "f");
  legend->AddEntry(all_m1_histo->GetName(), legendEntry3.c_str(), "f");
  legend->AddEntry(all_m2_histo->GetName(), legendEntry4.c_str(), "f");
  legend->AddEntry(all_m3_histo->GetName(), legendEntry5.c_str(), "f");
  legend->AddEntry(EFT_signal->GetName(), legendEntry6.c_str(), "l");
  legend->Draw("SAME");

  Canvas->Print(Name.c_str());


  double EFTmax = 3;
  double pas = EFTmax*2/EFT_histo_elec[0]->GetNbinsX();

  TH1D* chiCarre = new TH1D ("Chi Carre", "" , EFT_histo_elec[0]->GetNbinsX(), -EFTmax, EFTmax);





  for (int factor = 1 ; factor <= EFT_histo_elec[0]->GetNbinsX() ; factor++)
    {
      double normalisation_muon = 0;
      double normalisation_elec = 0;
      for (int p = 1; p <= nbins ; p++)
        {
          normalisation_elec += pow(ASIMOV_elec->GetBinContent(p) - (syst_lumi*(EFT_histo_elec[p-1]->GetBinContent(factor)*singletop_elec->GetBinContent(p)*syst_singletop + ttbar_tw_elec->GetBinContent(p)*syst_ttbar_tw + w_z_gamma_elec->GetBinContent(p)*syst_w_z_gamma)+ QCD_elec->GetBinContent(p)*syst_QCD),2)
        /(ASIMOV_elec->GetBinContent(p));
          normalisation_muon += pow(ASIMOV_muon->GetBinContent(p) - (syst_lumi*(EFT_histo_muon[p-1]->GetBinContent(factor)*singletop_muon->GetBinContent(p)*syst_singletop + ttbar_tw_muon->GetBinContent(p)*syst_ttbar_tw + w_z_gamma_muon->GetBinContent(p)*syst_w_z_gamma)+ QCD_muon->GetBinContent(p)*syst_QCD),2)
        /(ASIMOV_muon->GetBinContent(p));
        }
      chiCarre->SetBinContent(factor, (normalisation_elec+normalisation_muon)/nbins);
    }

double minimum_chicarre = chiCarre->GetMinimum();
for (int factor = 1 ; factor <= EFT_histo_elec[0]->GetNbinsX() ; factor++)
  {
    double current_value = chiCarre->GetBinContent(factor);
    chiCarre->SetBinContent(factor,current_value-minimum_chicarre);
  }

  chiCarre->SetXTitle("C_{tW}^{I}");
  chiCarre->SetYTitle("#\chi^{2}");

  chiCarre->Draw();
  Canvas->Print("chi_carre.pdf");

double bin_left = 99;
double bin_right = 99;


for (double factor = 1 ; factor < EFT_histo_elec[0]->GetNbinsX() ; factor++)
{
  double value_hold = chiCarre->GetBinContent(factor);
  double value_next = chiCarre->GetBinContent(factor+1);

  if(value_hold<1.0 && value_next>1.0 )
  {
    if (abs(50-((2*factor+1)/2)) < abs(50-bin_right))
      {
        bin_right = (2*factor+1)/2;
      }
    /*cout<<"value hold ="<<value_hold<<endl;
    cout<<"value_next ="<<value_next<<endl;*/
  }
  if(value_hold>1.0 && value_next<1.0)
  {
    if (abs(50-((2*factor+1)/2)) < abs(50-bin_left))
      {
        bin_left = (2*factor+1)/2;
      }



    /*cout<<"value hold ="<<value_hold<<endl;
    cout<<"value_next ="<<value_next<<endl;*/
  }

}
  cout<<"bin_left="<<bin_left<<endl;
  cout<<"bin_right="<<bin_right<<endl;

  double incertitude_left = -EFTmax+bin_left*pas;
  double incertitude_right = -EFTmax+bin_right*pas;
  double moyenne_incert = (abs(incertitude_left)+abs(incertitude_right))/2;
  cout<<"incertitude_left ="<<incertitude_left<<endl;
  cout<<"incertitude_right ="<<incertitude_right<<endl;
  cout<<"moyenne_incert ="<<moyenne_incert<<endl;




  return moyenne_incert;
}

void Compare_3Histos(TTree* t1, TTree* t2, TTree* t3, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string legendEntry3, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "Histo_2");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "Histo_3");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);

  double a = Histo_1->Integral();
  double b = Histo_2->Integral();
  double c = Histo_3->Integral();
  Histo_1->Scale(1/a);
  Histo_2->Scale(1/b);
  if (c>0) Histo_3->Scale(1/c);
  cout << "a="<<a<<" b="<<b<<" c="<<c<<endl;

  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();
  if (c>0) max = (max>Histo_3->GetMaximum()) ? max : Histo_3->GetMaximum();
  cout << "max="<<max<<endl;

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");

  //Canvas->SetLogx();
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
  if (legendX == "Top Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.0001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  if (legendX == "W Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

  if (c>0){
    Histo_3->SetLineColor(kOrange);
    Histo_3->SetLineWidth(2);
    Histo_3->Draw("SAME");
  }

 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.75;
	 lx1 = 0.5;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.75;
	 ly0 = 0.75;
	 lx1 = 0.99;
	 ly1 = 0.99;
  }

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");
  legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(), "l");
  if (c>0) legend->AddEntry(Histo_3->GetName(), legendEntry3.c_str(), "l");

  legend->Draw("SAME");

  Canvas->Print(Name.c_str());

  cout << "Histo1 mean: "<<Histo_1->GetMean()<<endl;
  cout << "Histo2 mean: "<<Histo_2->GetMean()<<endl;
  if (c>0) cout << "Histo3 mean: "<<Histo_3->GetMean()<<endl;

}

void Compare_1Histos(TTree* t1, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  Histo_1->SetStats(kFALSE);

  //double max = Histo_1->GetMaximum();
  double max = 3000;
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
  if (legendX =="electron isolation" || legendX =="muon isolation")
  {
    Histo_1->SetAxisRange(0.5,max*1.1,"Y");
    Canvas->SetLogy();
    Canvas->SetLogx();
  }

  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  double lx0, ly0, lx1, ly1;
   if (legendPlace=="legendUpLeft"){
 	 lx0 = 0.2;
 	 ly0 = 0.75;
 	 lx1 = 0.5;
 	 ly1 = 0.95;
   }
    if (legendPlace=="legendUpRight"){
 	 lx0 = 0.75;
 	 ly0 = 0.75;
 	 lx1 = 0.99;
 	 ly1 = 0.99;
   }

   TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
   legend->SetFillColor(kWhite);
   legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");

   legend->Draw("SAME");

   Canvas->Print(Name.c_str());
}

void Compare_2Histos(TTree* t1, TTree* t2, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "Histo_2");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);

  double a = Histo_1->Integral();
  double b = Histo_2->Integral();
  Histo_1->Scale(1/a);
  Histo_2->Scale(1/b);
  cout<<legendX<<endl;
  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
  if (legendX == "Top Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.0001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  if (legendX == "W Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.75;
	 lx1 = 0.5;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.75;
	 ly0 = 0.75;
	 lx1 = 0.99;
	 ly1 = 0.99;
  }

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");
  legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(), "l");
  legend->Draw("SAME");

  Canvas->Print(Name.c_str());

  cout << "Histo1 mean: "<<Histo_1->GetMean()<<endl;
  cout << "Histo2 mean: "<<Histo_2->GetMean()<<endl;
}

void Compare_4Histos(TTree* t1, TTree* t2, TTree* t3, TTree* t4, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string legendEntry3, string legendEntry4, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "Histo_2");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "Histo_3");
  TH1D* Histo_4 = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection, "Histo_4");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);
  Histo_4->SetStats(kFALSE);

  double a = Histo_1->Integral();
  double b = Histo_2->Integral();
  double c = Histo_3->Integral();
  double d = Histo_4->Integral();

  Histo_1->Scale(1/a);
  Histo_2->Scale(1/b);
  if (c>0) Histo_3->Scale(1/c);
  Histo_4->Scale(1/d);
  cout << "a="<<a<<" b="<<b<<" c="<<c<<endl;

  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();
  if (c>0) max = (max>Histo_3->GetMaximum()) ? max : Histo_3->GetMaximum();
  if (d>0) max = (max>Histo_4->GetMaximum()) ? max : Histo_4->GetMaximum();
  cout << "max="<<max<<endl;

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  //Canvas->SetLogx();
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
  if (legendX == "Top Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.0001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  if (legendX == "W Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

  if (c>0){
    Histo_3->SetLineColor(kOrange);
    Histo_3->SetLineWidth(2);
    Histo_3->Draw("SAME");
  }

  Histo_4->SetLineColor(kGreen);
  Histo_4->SetLineWidth(2);
  Histo_4->Draw("SAME");

 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.75;
	 lx1 = 0.5;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.75;
	 ly0 = 0.75;
	 lx1 = 0.99;
	 ly1 = 0.99;
  }

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");
  legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(), "l");
  legend->AddEntry(Histo_3->GetName(), legendEntry3.c_str(), "l");
  legend->AddEntry(Histo_4->GetName(), legendEntry4.c_str(), "l");
  legend->Draw("SAME");

  Canvas->Print(Name.c_str());

  cout << "Histo1 mean: "<<Histo_1->GetMean()<<endl;
  cout << "Histo2 mean: "<<Histo_2->GetMean()<<endl;
  if (c>0) cout << "Histo3 mean: "<<Histo_3->GetMean()<<endl;

}
void Stacked_histo(TTree* t1, TTree* t2, TTree* t3, TTree* t4 , TTree* t5, TTree* t6, TTree* t7, TTree* t8, TTree* t9, TTree* t10, TTree* t11, TTree* t12, TTree* t13, TTree* t14, TTree* t15, TTree* t16, TTree* t19, TTree* t20, TTree* t21, TTree* t22, TTree* t23, TTree* t24,
   TTree* t25, TTree* t26, TTree* t27, TTree* t28, TTree* t29,
   string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string Name ,double param0, double param1, string region){

  double lumi = 41500;
  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "Histo_2");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "Histo_3");
  TH1D* Histo_4 = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection, "Histo_4");
  TH1D* Histo_5 = GetHistoWeight(t5, variable, nbins, xmin, xmax, selection, "Histo_5");
  TH1D* Histo_6 = GetHistoWeight(t6, variable, nbins, xmin, xmax, selection, "Histo_6");
  TH1D* Histo_7 = GetHistoWeight(t7, variable, nbins, xmin, xmax, selection, "Histo_7");
  TH1D* Histo_8 = GetHistoWeight(t8, variable, nbins, xmin, xmax, selection, "Histo_8");
  TH1D* Histo_9 = GetHistoWeight(t9, variable, nbins, xmin, xmax, selection, "Histo_9");
  TH1D* Histo_10 = GetHistoWeight(t10, variable, nbins, xmin, xmax, selection, "Histo_10");
  TH1D* Histo_11 = GetHistoWeight(t11, variable, nbins, xmin, xmax, selection, "Histo_11");
  TH1D* Histo_12 = GetHistoWeight(t12, variable, nbins, xmin, xmax, selection, "Histo_12");
  TH1D* Histo_13 = GetHistoWeight(t13, variable, nbins, xmin, xmax, selection, "Histo_13");
  TH1D* Histo_14 = GetHistoWeight(t14, variable, nbins, xmin, xmax, selection, "Histo_14");
  TH1D* Histo_15 = GetHistoWeight(t15, variable, nbins, xmin, xmax, selection, "Histo_15");
  TH1D* Histo_16 = GetHistoWeight(t16, variable, nbins, xmin, xmax, selection, "Histo_16");
  /*TH1D* Histo_17 = GetHistoWeight(t17, variable, nbins, xmin, xmax, selection, "Histo_17");
  TH1D* Histo_18 = GetHistoWeight(t18, variable, nbins, xmin, xmax, selection, "Histo_18");*/
  TH1D* Histo_19 = GetHistoWeight(t19, variable, nbins, xmin, xmax, selection, "Histo_19");
  TH1D* Histo_20 = GetHistoWeight(t20, variable, nbins, xmin, xmax, selection, "Histo_20");
  TH1D* Histo_21 = GetHistoWeight(t21, variable, nbins, xmin, xmax, selection, "Histo_21");
  TH1D* Histo_22 = GetHistoWeight(t22, variable, nbins, xmin, xmax, selection, "Histo_22");
  TH1D* Histo_23 = GetHistoWeight(t23, variable, nbins, xmin, xmax, selection, "Histo_23");
  TH1D* Histo_24 = GetHistoWeight(t24, variable, nbins, xmin, xmax, selection, "Histo_24");
  TH1D* Histo_25 = GetHistoWeight(t25, variable, nbins, xmin, xmax, selection, "Histo_25");
  TH1D* Histo_26 = GetHistoWeight(t26, variable, nbins, xmin, xmax, selection, "Histo_26");
  TH1D* Histo_27 = GetHistoWeight(t27, variable, nbins, xmin, xmax, selection, "Histo_27");
  TH1D* Histo_28 = GetHistoWeight(t28, variable, nbins, xmin, xmax, selection, "Histo_28");
  TH1D* DATA = GetHistoWeight(t29, variable, nbins, xmin, xmax, selection, "Histo_29");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);
  Histo_4->SetStats(kFALSE);
  Histo_5->SetStats(kFALSE);
  Histo_6->SetStats(kFALSE);
  Histo_7->SetStats(kFALSE);
  Histo_8->SetStats(kFALSE);
  Histo_9->SetStats(kFALSE);
  Histo_10->SetStats(kFALSE);
  Histo_11->SetStats(kFALSE);
  Histo_12->SetStats(kFALSE);
  Histo_13->SetStats(kFALSE);
  Histo_14->SetStats(kFALSE);
  Histo_15->SetStats(kFALSE);
  Histo_16->SetStats(kFALSE);
  /*Histo_17->SetStats(kFALSE);
  Histo_18->SetStats(kFALSE);*/
  Histo_19->SetStats(kFALSE);
  Histo_20->SetStats(kFALSE);
  Histo_21->SetStats(kFALSE);
  Histo_22->SetStats(kFALSE);
  Histo_23->SetStats(kFALSE);
  Histo_24->SetStats(kFALSE);
  Histo_25->SetStats(kFALSE);
  Histo_26->SetStats(kFALSE);
  Histo_27->SetStats(kFALSE);
  Histo_28->SetStats(kFALSE);
  DATA->SetStats(kFALSE);


  double events_MC_t_top = 5982064;
  double events_MC_t_antitop = 3675910;
  double events_MC_DY_50 = 491250940213;
  double events_MC_DY_1050 = 315866;
  double events_MC_semilep = 32414370098;
  double events_MC_tW_antitop = 270762750.172;
  double events_MC_tW_top = 277241050.839;
  double events_MC_W_jets = 29981320;
  double events_MC_dilep = 4984888995.13;
  double events_MC_QCD_EM_20to30 = 11212816.767;
  double events_MC_QCD_EM_30to50 = 14766010 ;
  double events_MC_QCD_EM_50to80 = 10288889 ;
  double events_MC_QCD_EM_80to120 = 8731568.95822;
  double events_MC_QCD_EM_120to170 = 8351696.45656 ;
  double events_MC_QCD_EM_170to300 = 3433364 ;
  double events_MC_QCD_EM_300toInf = 2874295;
  double events_MC_QCD_Mu_15to20 = 5860867.35455;
  double events_MC_QCD_Mu_20to30 = 27812825.2452;
  double events_MC_QCD_Mu_30to50 = 29030324 ;
  double events_MC_QCD_Mu_50to80 = 24068613 ;
  double events_MC_QCD_Mu_80to120 = 23251833.189;
  double events_MC_QCD_Mu_120to170 = 20775373.413 ;
  double events_MC_QCD_Mu_170to300 = 46170668;
  double events_MC_QCD_Mu_300to470 = 17336667.2501;
  double events_MC_QCD_Mu_470to600 = 23515446.0718;
  double events_MC_QCD_Mu_600to800 = 17263676.523 ;
  double events_MC_QCD_Mu_800to1000 = 16520155;
  double events_MC_QCD_Mu_1000toInf = 11148422;

  double cross_MC_t_top = 136.02;
  double cross_MC_t_antitop = 80.95;
  double cross_MC_DY_50 = 6025.2;
  double cross_MC_DY_1050 = 22635.1;
  double cross_MC_semilep = 365.3;
  double cross_MC_tW_antitop = 35.9;
  double cross_MC_tW_top = 35.9;
  double cross_MC_W_jets = 61526.7;
  double cross_MC_dilep = 88.2;
  double cross_MC_QCD_EM_20to30 = 4948840 ;
  double cross_MC_QCD_EM_30to50 = 6324800 ;
  double cross_MC_QCD_EM_50to80 = 1806336;
  double cross_MC_QCD_EM_80to120 = 380538;
  double cross_MC_QCD_EM_120to170 = 66634.28 ;
  double cross_MC_QCD_EM_170to300 = 20859 ;
  double cross_MC_QCD_EM_300toInf = 1350 ;
  double cross_MC_QCD_Mu_15to20 = 3336011;
  double cross_MC_QCD_Mu_20to30 = 3987854.9 ;
  double cross_MC_QCD_Mu_30to50 = 1705381 ;
  double cross_MC_QCD_Mu_50to80 = 395178 ;
  double cross_MC_QCD_Mu_80to120 = 106889.4;
  double cross_MC_QCD_Mu_120to170 = 23773.61 ;
  double cross_MC_QCD_Mu_170to300 = 8292.982 ;
  double cross_MC_QCD_Mu_300to470 = 797.35269;
  double cross_MC_QCD_Mu_470to600 = 56.588336 ;
  double cross_MC_QCD_Mu_600to800 = 25.09505908 ;
  double cross_MC_QCD_Mu_800to1000 = 4.707368272;
  double cross_MC_QCD_Mu_1000toInf = 1.62131692;

  Histo_1->Scale(1.0/events_MC_t_top*cross_MC_t_top*lumi);
  Histo_2->Scale(1.0/events_MC_t_antitop*cross_MC_t_antitop*lumi);
  Histo_3->Scale(1.0/events_MC_DY_50*cross_MC_DY_50*lumi);
  Histo_4->Scale(1.0/events_MC_DY_1050*cross_MC_DY_1050*lumi);
  Histo_5->Scale(1.0/events_MC_semilep*cross_MC_semilep*lumi);
  Histo_6->Scale(1.0/events_MC_tW_antitop*cross_MC_tW_antitop*lumi);
  Histo_7->Scale(1.0/events_MC_tW_top*cross_MC_tW_top*lumi);
  Histo_8->Scale(1/events_MC_W_jets*cross_MC_W_jets*lumi);
  Histo_9->Scale(1.0/events_MC_dilep*cross_MC_dilep*lumi);
  Histo_10->Scale(1.0/events_MC_QCD_EM_20to30*cross_MC_QCD_EM_20to30*lumi);
  Histo_11->Scale(1.0/events_MC_QCD_EM_30to50*cross_MC_QCD_EM_30to50*lumi);
  Histo_12->Scale(1.0/events_MC_QCD_EM_50to80*cross_MC_QCD_EM_50to80*lumi);
  Histo_13->Scale(1.0/events_MC_QCD_EM_80to120*cross_MC_QCD_EM_80to120*lumi);
  Histo_14->Scale(1.0/events_MC_QCD_EM_120to170*cross_MC_QCD_EM_120to170*lumi);
  Histo_15->Scale(1.0/events_MC_QCD_EM_170to300*cross_MC_QCD_EM_170to300*lumi);
  Histo_16->Scale(1.0/events_MC_QCD_EM_300toInf*cross_MC_QCD_EM_300toInf*lumi);
  /*Histo_17->Scale(1.0/events_MC_QCD_Mu_15to20*cross_MC_QCD_Mu_15to20*lumi);
  Histo_18->Scale(1.0/events_MC_QCD_Mu_20to30*cross_MC_QCD_Mu_20to30*lumi);*/
  Histo_19->Scale(1.0/events_MC_QCD_Mu_30to50*cross_MC_QCD_Mu_30to50*lumi);
  Histo_20->Scale(1.0/events_MC_QCD_Mu_50to80*cross_MC_QCD_Mu_50to80*lumi);
  Histo_21->Scale(1.0/events_MC_QCD_Mu_80to120*cross_MC_QCD_Mu_80to120*lumi);
  Histo_22->Scale(1.0/events_MC_QCD_Mu_120to170*cross_MC_QCD_Mu_120to170*lumi);
  Histo_23->Scale(1.0/events_MC_QCD_Mu_170to300*cross_MC_QCD_Mu_170to300*lumi);
  Histo_24->Scale(1.0/events_MC_QCD_Mu_300to470*cross_MC_QCD_Mu_300to470*lumi);
  Histo_25->Scale(1.0/events_MC_QCD_Mu_470to600*cross_MC_QCD_Mu_470to600*lumi);
  Histo_26->Scale(1.0/events_MC_QCD_Mu_600to800*cross_MC_QCD_Mu_600to800*lumi);
  Histo_27->Scale(1.0/events_MC_QCD_Mu_800to1000*cross_MC_QCD_Mu_800to1000*lumi);
  Histo_28->Scale(1.0/events_MC_QCD_Mu_1000toInf*cross_MC_QCD_Mu_1000toInf*lumi);

  //TH1D *test = new TH1D ("stack single top", "", nbins, xmin, xmax );
  TH1D *singletop = new TH1D ("stack single top", "", nbins, xmin, xmax );
  TH1D *ttbar_tw= new TH1D ("stack ttbar_tw", "", nbins, xmin, xmax );
  TH1D *QCD = new TH1D ("stack QCD", "", nbins, xmin, xmax );
  TH1D *w_z_gamma = new TH1D ("stack w_z_gamma", "", nbins, xmin, xmax );

  singletop->Add(Histo_1);
  singletop->Add(Histo_2);

  ttbar_tw->Add(Histo_5);
  ttbar_tw->Add(Histo_9);
  ttbar_tw->Add(Histo_6);
  ttbar_tw->Add(Histo_7);


  w_z_gamma->Add(Histo_3);
  w_z_gamma->Add(Histo_4);
  w_z_gamma->Add(Histo_8);

  QCD->Add(Histo_10);
  QCD->Add(Histo_11);
  QCD->Add(Histo_12);
  QCD->Add(Histo_13);
  QCD->Add(Histo_14);
  QCD->Add(Histo_15);
  QCD->Add(Histo_16);
  /*QCD->Add(Histo_17);
  QCD->Add(Histo_18);*/
  QCD->Add(Histo_19);
  QCD->Add(Histo_20);
  QCD->Add(Histo_21);
  QCD->Add(Histo_22);
  QCD->Add(Histo_23);
  QCD->Add(Histo_24);
  QCD->Add(Histo_25);
  QCD->Add(Histo_26);
  QCD->Add(Histo_27);
  QCD->Add(Histo_28);

  TH1D *all_others = new TH1D ("all_others","", nbins, xmin, xmax);

  all_others->Add(singletop);
  all_others->Add(ttbar_tw);
  all_others->Add(w_z_gamma);

  double integ_all_others = all_others->Integral(1,20);
  double integ_QCD = QCD->Integral(1,20);
  double integ_data = DATA->Integral(1,20);

  double constrain1 = integ_all_others/integ_data;

  cout<<"integral QCD ="<<integ_QCD<<endl;

  cout<<"all others proportion = "<<integ_all_others/integ_data<<endl;
  cout<<"integral DATA ="<<integ_data<<endl;

  /*double newscale0 = integ_data*param0/integ_all_others;
  double newscale1 = integ_data*param1/integ_QCD;
  cout<<"newscale0 = "<< newscale0<<endl;
  cout<<"newscale1 = "<< newscale1<<endl;
  singletop->Scale(newscale0);
  ttbar_tw->Scale(newscale0);
  w_z_gamma->Scale(newscale0);
  QCD->Scale(newscale1);*/

  TH1D *all_histo = new TH1D ("multijet", "multijet", nbins, xmin, xmax );
  TH1D *all_m1_histo = new TH1D ("w_z_gamma", "w_z_gamma", nbins, xmin, xmax );
  TH1D *all_m2_histo = new TH1D ("ttbar_tw", "ttbar_tw", nbins, xmin, xmax );
  TH1D *all_m3_histo = new TH1D ("singletop", "singletop", nbins, xmin, xmax );


  all_histo->Add(singletop);
  all_histo->Add(ttbar_tw);
  all_histo->Add(w_z_gamma);

  all_histo->Add(QCD);


  all_m1_histo->Add(singletop);
  all_m1_histo->Add(ttbar_tw);
  all_m1_histo->Add(w_z_gamma);

  all_m2_histo->Add(singletop);
  all_m2_histo->Add(ttbar_tw);

  all_m3_histo->Add(singletop);

  all_histo->SetLineColor(kBlack);
  all_histo->SetFillColor(kGray);

  all_m1_histo->SetLineColor(kBlack);
  all_m1_histo->SetFillColor(8);

  all_m2_histo->SetLineColor(kBlack);
  all_m2_histo->SetFillColor(kOrange);

  all_m3_histo->SetLineColor(kBlack);
  all_m3_histo->SetFillColor(kRed);

  DATA->SetLineWidth(2);
  DATA->SetLineColor(kBlack);
  all_histo->SetXTitle(legendX.c_str());
  all_histo->SetYTitle(legendY.c_str());

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  double max = all_histo->GetMaximum();
  all_histo->SetAxisRange(0,max*1.2,"Y");
  //Canvas->SetLogx();

  //Canvas->SetLogy();
  all_histo->Draw("HIST");
  all_m1_histo->Draw("SAME HIST");
  all_m2_histo->Draw("SAME HIST");
  all_m3_histo->Draw("SAME HIST");
  DATA->Draw("E SAME");


  string legendtitle1 = region +" Region";
  float lx0 = 0.65;
  float ly0 = 0.65;
  float lx1 = 0.99;
  float ly1 = 0.99;

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle1.c_str());
  string legendEntry1 = "DATA";
  string legendEntry2 = "Multijet";
  string legendEntry3 = "W/Z/#gamma* + jets";
  string legendEntry4 = "t#bar{t}/tw";
  string legendEntry5 = "Single top t-channel";




  legend->SetFillColor(kWhite);
  legend->AddEntry(DATA->GetName(), legendEntry1.c_str(), "l");
  legend->AddEntry(all_histo->GetName(), legendEntry2.c_str(), "f");
  legend->AddEntry(all_m1_histo->GetName(), legendEntry3.c_str(), "f");
  legend->AddEntry(all_m2_histo->GetName(), legendEntry4.c_str(), "f");
  legend->AddEntry(all_m3_histo->GetName(), legendEntry5.c_str(), "f");
    legend->Draw("SAME");

  Canvas->Print(Name.c_str());
}

void Stacked_histo_Fit(TTree* t1, TTree* t2, TTree* t3, TTree* t4 , TTree* t5, TTree* t6, TTree* t7, TTree* t8, TTree* t9, TTree* t10, TTree* t11, TTree* t12, TTree* t13, TTree* t14, TTree* t15, TTree* t16, TTree* t19, TTree* t20, TTree* t21, TTree* t22, TTree* t23, TTree* t24,
   TTree* t25, TTree* t26, TTree* t27, TTree* t28, TTree* t29,
   string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string Name){

  double lumi = 41500;
  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "Histo_2");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "Histo_3");
  TH1D* Histo_4 = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection, "Histo_4");
  TH1D* Histo_5 = GetHistoWeight(t5, variable, nbins, xmin, xmax, selection, "Histo_5");
  TH1D* Histo_6 = GetHistoWeight(t6, variable, nbins, xmin, xmax, selection, "Histo_6");
  TH1D* Histo_7 = GetHistoWeight(t7, variable, nbins, xmin, xmax, selection, "Histo_7");
  TH1D* Histo_8 = GetHistoWeight(t8, variable, nbins, xmin, xmax, selection, "Histo_8");
  TH1D* Histo_9 = GetHistoWeight(t9, variable, nbins, xmin, xmax, selection, "Histo_9");
  TH1D* Histo_10 = GetHistoWeight(t10, variable, nbins, xmin, xmax, selection, "Histo_10");
  TH1D* Histo_11 = GetHistoWeight(t11, variable, nbins, xmin, xmax, selection, "Histo_11");
  TH1D* Histo_12 = GetHistoWeight(t12, variable, nbins, xmin, xmax, selection, "Histo_12");
  TH1D* Histo_13 = GetHistoWeight(t13, variable, nbins, xmin, xmax, selection, "Histo_13");
  TH1D* Histo_14 = GetHistoWeight(t14, variable, nbins, xmin, xmax, selection, "Histo_14");
  TH1D* Histo_15 = GetHistoWeight(t15, variable, nbins, xmin, xmax, selection, "Histo_15");
  TH1D* Histo_16 = GetHistoWeight(t16, variable, nbins, xmin, xmax, selection, "Histo_16");
  //TH1D* Histo_17 = GetHistoWeight(t17, variable, nbins, xmin, xmax, selection, "Histo_17");
  //TH1D* Histo_18 = GetHistoWeight(t18, variable, nbins, xmin, xmax, selection, "Histo_18");
  TH1D* Histo_19 = GetHistoWeight(t19, variable, nbins, xmin, xmax, selection, "Histo_19");
  TH1D* Histo_20 = GetHistoWeight(t20, variable, nbins, xmin, xmax, selection, "Histo_20");
  TH1D* Histo_21 = GetHistoWeight(t21, variable, nbins, xmin, xmax, selection, "Histo_21");
  TH1D* Histo_22 = GetHistoWeight(t22, variable, nbins, xmin, xmax, selection, "Histo_22");
  TH1D* Histo_23 = GetHistoWeight(t23, variable, nbins, xmin, xmax, selection, "Histo_23");
  TH1D* Histo_24 = GetHistoWeight(t24, variable, nbins, xmin, xmax, selection, "Histo_24");
  TH1D* Histo_25 = GetHistoWeight(t25, variable, nbins, xmin, xmax, selection, "Histo_25");
  TH1D* Histo_26 = GetHistoWeight(t26, variable, nbins, xmin, xmax, selection, "Histo_26");
  TH1D* Histo_27 = GetHistoWeight(t27, variable, nbins, xmin, xmax, selection, "Histo_27");
  TH1D* Histo_28 = GetHistoWeight(t28, variable, nbins, xmin, xmax, selection, "Histo_28");
  TH1D* DATA = GetHistoWeight(t29, variable, nbins, xmin, xmax, selection, "Histo_29");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);
  Histo_4->SetStats(kFALSE);
  Histo_5->SetStats(kFALSE);
  Histo_6->SetStats(kFALSE);
  Histo_7->SetStats(kFALSE);
  Histo_8->SetStats(kFALSE);
  Histo_9->SetStats(kFALSE);
  Histo_10->SetStats(kFALSE);
  Histo_11->SetStats(kFALSE);
  Histo_12->SetStats(kFALSE);
  Histo_13->SetStats(kFALSE);
  Histo_14->SetStats(kFALSE);
  Histo_15->SetStats(kFALSE);
  Histo_16->SetStats(kFALSE);
  //Histo_17->SetStats(kFALSE);
  //Histo_18->SetStats(kFALSE);
  Histo_19->SetStats(kFALSE);
  Histo_20->SetStats(kFALSE);
  Histo_21->SetStats(kFALSE);
  Histo_22->SetStats(kFALSE);
  Histo_23->SetStats(kFALSE);
  Histo_24->SetStats(kFALSE);
  Histo_25->SetStats(kFALSE);
  Histo_26->SetStats(kFALSE);
  Histo_27->SetStats(kFALSE);
  Histo_28->SetStats(kFALSE);
  DATA->SetStats(kFALSE);

  double events_MC_t_top = 5982064;
  double events_MC_t_antitop = 3675910;
  double events_MC_DY_50 = 491250940213;
  double events_MC_DY_1050 = 315866;
  double events_MC_semilep = 32414370098;
  double events_MC_tW_antitop = 270762750.172;
  double events_MC_tW_top = 277241050.839;
  double events_MC_W_jets = 29981320;
  double events_MC_dilep = 4984888995.13;
  double events_MC_QCD_EM_20to30 = 11212816.767;
  double events_MC_QCD_EM_30to50 = 14766010 ;
  double events_MC_QCD_EM_50to80 = 10288889 ;
  double events_MC_QCD_EM_80to120 = 8731568.95822;
  double events_MC_QCD_EM_120to170 = 8351696.45656 ;
  double events_MC_QCD_EM_170to300 = 3433364 ;
  double events_MC_QCD_EM_300toInf = 2874295;
  //double events_MC_QCD_Mu_15to20 = 5860867.35455;
  //double events_MC_QCD_Mu_20to30 = 27812825.2452;
  double events_MC_QCD_Mu_30to50 = 29030324 ;
  double events_MC_QCD_Mu_50to80 = 24068613 ;
  double events_MC_QCD_Mu_80to120 = 23251833.189;
  double events_MC_QCD_Mu_120to170 = 20775373.413 ;
  double events_MC_QCD_Mu_170to300 = 46170668;
  double events_MC_QCD_Mu_300to470 = 17336667.2501;
  double events_MC_QCD_Mu_470to600 = 23515446.0718;
  double events_MC_QCD_Mu_600to800 = 17263676.523 ;
  double events_MC_QCD_Mu_800to1000 = 16520155;
  double events_MC_QCD_Mu_1000toInf = 11148422;

  double cross_MC_t_top = 136.02;
  double cross_MC_t_antitop = 80.95;
  double cross_MC_DY_50 = 6025.2;
  double cross_MC_DY_1050 = 22635.1;
  double cross_MC_semilep = 365.3;
  double cross_MC_tW_antitop = 35.9;
  double cross_MC_tW_top = 35.9;
  double cross_MC_W_jets = 61526.7;
  double cross_MC_dilep = 88.2;
  double cross_MC_QCD_EM_20to30 = 4948840 ;
  double cross_MC_QCD_EM_30to50 = 6324800 ;
  double cross_MC_QCD_EM_50to80 = 1806336;
  double cross_MC_QCD_EM_80to120 = 380538;
  double cross_MC_QCD_EM_120to170 = 66634.28 ;
  double cross_MC_QCD_EM_170to300 = 20859 ;
  double cross_MC_QCD_EM_300toInf = 1350 ;
  //double cross_MC_QCD_Mu_15to20 = 3336011;
  //double cross_MC_QCD_Mu_20to30 = 3987854.9 ;
  double cross_MC_QCD_Mu_30to50 = 1705381 ;
  double cross_MC_QCD_Mu_50to80 = 395178 ;
  double cross_MC_QCD_Mu_80to120 = 106889.4;
  double cross_MC_QCD_Mu_120to170 = 23773.61 ;
  double cross_MC_QCD_Mu_170to300 = 8292.982 ;
  double cross_MC_QCD_Mu_300to470 = 797.35269;
  double cross_MC_QCD_Mu_470to600 = 56.588336 ;
  double cross_MC_QCD_Mu_600to800 = 25.09505908 ;
  double cross_MC_QCD_Mu_800to1000 = 4.707368272;
  double cross_MC_QCD_Mu_1000toInf = 1.62131692;

  Histo_1->Scale(1.0/events_MC_t_top*cross_MC_t_top*lumi);
  Histo_2->Scale(1.0/events_MC_t_antitop*cross_MC_t_antitop*lumi);
  Histo_3->Scale(1.0/events_MC_DY_50*cross_MC_DY_50*lumi);
  Histo_4->Scale(1.0/events_MC_DY_1050*cross_MC_DY_1050*lumi);
  Histo_5->Scale(1.0/events_MC_semilep*cross_MC_semilep*lumi);
  Histo_6->Scale(1.0/events_MC_tW_antitop*cross_MC_tW_antitop*lumi);
  Histo_7->Scale(1.0/events_MC_tW_top*cross_MC_tW_top*lumi);
  Histo_8->Scale(1/events_MC_W_jets*cross_MC_W_jets*lumi);
  Histo_9->Scale(1.0/events_MC_dilep*cross_MC_dilep*lumi);
  Histo_10->Scale(1.0/events_MC_QCD_EM_20to30*cross_MC_QCD_EM_20to30*lumi);
  Histo_11->Scale(1.0/events_MC_QCD_EM_30to50*cross_MC_QCD_EM_30to50*lumi);
  Histo_12->Scale(1.0/events_MC_QCD_EM_50to80*cross_MC_QCD_EM_50to80*lumi);
  Histo_13->Scale(1.0/events_MC_QCD_EM_80to120*cross_MC_QCD_EM_80to120*lumi);
  Histo_14->Scale(1.0/events_MC_QCD_EM_120to170*cross_MC_QCD_EM_120to170*lumi);
  Histo_15->Scale(1.0/events_MC_QCD_EM_170to300*cross_MC_QCD_EM_170to300*lumi);
  Histo_16->Scale(1.0/events_MC_QCD_EM_300toInf*cross_MC_QCD_EM_300toInf*lumi);
  //Histo_17->Scale(1.0/events_MC_QCD_Mu_15to20*cross_MC_QCD_Mu_15to20*lumi);
  //Histo_18->Scale(1.0/events_MC_QCD_Mu_20to30*cross_MC_QCD_Mu_20to30*lumi);
  Histo_19->Scale(1.0/events_MC_QCD_Mu_30to50*cross_MC_QCD_Mu_30to50*lumi);
  Histo_20->Scale(1.0/events_MC_QCD_Mu_50to80*cross_MC_QCD_Mu_50to80*lumi);
  Histo_21->Scale(1.0/events_MC_QCD_Mu_80to120*cross_MC_QCD_Mu_80to120*lumi);
  Histo_22->Scale(1.0/events_MC_QCD_Mu_120to170*cross_MC_QCD_Mu_120to170*lumi);
  Histo_23->Scale(1.0/events_MC_QCD_Mu_170to300*cross_MC_QCD_Mu_170to300*lumi);
  Histo_24->Scale(1.0/events_MC_QCD_Mu_300to470*cross_MC_QCD_Mu_300to470*lumi);
  Histo_25->Scale(1.0/events_MC_QCD_Mu_470to600*cross_MC_QCD_Mu_470to600*lumi);
  Histo_26->Scale(1.0/events_MC_QCD_Mu_600to800*cross_MC_QCD_Mu_600to800*lumi);
  Histo_27->Scale(1.0/events_MC_QCD_Mu_800to1000*cross_MC_QCD_Mu_800to1000*lumi);
  Histo_28->Scale(1.0/events_MC_QCD_Mu_1000toInf*cross_MC_QCD_Mu_1000toInf*lumi);

  //TH1D *test = new TH1D ("stack single top", "", nbins, xmin, xmax );
  TH1D *singletop = new TH1D ("stack single top", "", nbins, xmin, xmax );
  TH1D *ttbar_tw= new TH1D ("stack ttbar_tw", "", nbins, xmin, xmax );
  TH1D *QCD = new TH1D ("stack QCD", "", nbins, xmin, xmax );
  TH1D *w_z_gamma = new TH1D ("stack w_z_gamma", "", nbins, xmin, xmax );

  singletop->Add(Histo_1);
  singletop->Add(Histo_2);

  ttbar_tw->Add(Histo_5);
  ttbar_tw->Add(Histo_9);
  ttbar_tw->Add(Histo_6);
  ttbar_tw->Add(Histo_7);


  w_z_gamma->Add(Histo_3);
  w_z_gamma->Add(Histo_4);
  w_z_gamma->Add(Histo_8);

  QCD->Add(Histo_10);
  QCD->Add(Histo_11);
  QCD->Add(Histo_12);
  QCD->Add(Histo_13);
  QCD->Add(Histo_14);
  QCD->Add(Histo_15);
  QCD->Add(Histo_16);
  //QCD->Add(Histo_17);
  //QCD->Add(Histo_18);
  QCD->Add(Histo_19);
  QCD->Add(Histo_20);
  QCD->Add(Histo_21);
  QCD->Add(Histo_22);
  QCD->Add(Histo_23);
  QCD->Add(Histo_24);
  QCD->Add(Histo_25);
  QCD->Add(Histo_26);
  QCD->Add(Histo_27);
  QCD->Add(Histo_28);

  TH1D *all_others = new TH1D ("all_others","", nbins, xmin, xmax);

  all_others->Add(singletop);
  all_others->Add(ttbar_tw);
  all_others->Add(w_z_gamma);

  double integ_all_others = all_others->Integral(1,20);
  double integ_QCD = QCD->Integral(1,20);
  double integ_data = DATA->Integral(1,20);

  double constrain1 = integ_all_others/integ_data;

  cout<<"integral QCD ="<<integ_QCD<<endl;

  cout<<"all others proportion = "<<integ_all_others/integ_data<<endl;
  cout<<"integral DATA ="<<integ_data<<endl;

  TObjArray *mc = new TObjArray(4);

  mc->Add(all_others);
  mc->Add(QCD);

  TFractionFitter* fit = new TFractionFitter(DATA,mc);
  fit->Constrain(0,constrain1-0.01,constrain1+0.01);
  fit->Constrain(1,0.0,0.5);

  fit->SetRangeX(1,20);
  Int_t status = fit->Fit();
  cout<<"fit status:"<<status<<endl;

  double param0, param1, param2, param3;
  double param0_error, param1_error, param2_error, param3_error;

  fit->GetResult(0,param0,param0_error);
  fit->GetResult(1,param1,param1_error);

  double newscale0 = integ_data*param0/integ_all_others;
  double newscale1 = integ_data*param1/integ_QCD;

  cout<<"newscale0="<<newscale0<<endl;
  cout<<"newscale1="<<newscale1<<endl;
  singletop->Scale(newscale0);
  ttbar_tw->Scale(newscale0);
  w_z_gamma->Scale(newscale0);
  QCD->Scale(newscale1);

  TH1D *all_histo = new TH1D ("multijet", "multijet", nbins, xmin, xmax );
  TH1D *all_m1_histo = new TH1D ("w_z_gamma", "w_z_gamma", nbins, xmin, xmax );
  TH1D *all_m2_histo = new TH1D ("ttbar_tw", "ttbar_tw", nbins, xmin, xmax );
  TH1D *all_m3_histo = new TH1D ("singletop", "singletop", nbins, xmin, xmax );


  all_histo->Add(singletop);
  all_histo->Add(ttbar_tw);
  all_histo->Add(w_z_gamma);
  all_histo->Add(QCD);

  all_m1_histo->Add(singletop);
  all_m1_histo->Add(ttbar_tw);
  all_m1_histo->Add(w_z_gamma);

  all_m2_histo->Add(singletop);
  all_m2_histo->Add(ttbar_tw);

  all_m3_histo->Add(singletop);

  all_histo->SetLineColor(kBlack);
  all_histo->SetFillColor(kGray);

  all_m1_histo->SetLineColor(kBlack);
  all_m1_histo->SetFillColor(8);

  all_m2_histo->SetLineColor(kBlack);
  all_m2_histo->SetFillColor(kOrange);

  all_m3_histo->SetLineColor(kBlack);
  all_m3_histo->SetFillColor(kRed);

  DATA->SetLineWidth(2);
  DATA->SetLineColor(kBlack);
  all_histo->SetXTitle(legendX.c_str());
  all_histo->SetYTitle(legendY.c_str());

  double max = all_histo->GetMaximum();
  all_histo->SetAxisRange(0,max*1.7,"Y");

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");

  all_histo->Draw("HIST");
  all_m1_histo->Draw("SAME HIST");
  all_m2_histo->Draw("SAME HIST");
  all_m3_histo->Draw("SAME HIST");
  DATA->Draw("E SAME");


  string legendtitle1 = "Signal Region";
  float lx0 = 0.65;
  float ly0 = 0.65;
  float lx1 = 0.99;
  float ly1 = 0.99;

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle1.c_str());
  string legendEntry1 = "DATA";
  string legendEntry2 = "Multijet";
  string legendEntry3 = "W/Z/#gamma* + jets";
  string legendEntry4 = "t#bar{t}/tw";
  string legendEntry5 = "Single top t-channel";



  legend->SetFillColor(kWhite);
  legend->AddEntry(DATA->GetName(), legendEntry1.c_str(), "l");
  legend->AddEntry(all_histo->GetName(), legendEntry2.c_str(), "f");
  legend->AddEntry(all_m1_histo->GetName(), legendEntry3.c_str(), "f");
  legend->AddEntry(all_m2_histo->GetName(), legendEntry4.c_str(), "f");
  legend->AddEntry(all_m3_histo->GetName(), legendEntry5.c_str(), "f");
  legend->Draw("SAME");

  Canvas->Print(Name.c_str());
}

void Stacked_histo_reversed(TTree* t1, TTree* t2, TTree* t3, TTree* t4 , TTree* t5, TTree* t6, TTree* t7, TTree* t8, TTree* t9, TTree* t10, TTree* t11,
   string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string Name, double param0, double param1, string lepton, string region){

  double lumi = 41500;

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_4 = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_5 = GetHistoWeight(t5, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_6 = GetHistoWeight(t6, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_7 = GetHistoWeight(t7, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_8 = GetHistoWeight(t8, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_9 = GetHistoWeight(t9, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_10 = GetHistoWeight(t10, variable, nbins, xmin, xmax, selection, "");
  TH1D* DATA = GetHistoWeight(t11, variable, nbins, xmin, xmax, selection, "");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);
  Histo_4->SetStats(kFALSE);
  Histo_5->SetStats(kFALSE);
  Histo_6->SetStats(kFALSE);
  Histo_7->SetStats(kFALSE);
  Histo_8->SetStats(kFALSE);
  Histo_9->SetStats(kFALSE);
  Histo_10->SetStats(kFALSE);

  DATA->SetStats(kFALSE);

  double events_MC_t_top = 5982064;
  double events_MC_t_antitop = 3675910;
  double events_MC_DY_50 = 491250940213;
  double events_MC_DY_1050 = 315866;
  double events_MC_semilep = 32414370098;
  double events_MC_tW_antitop = 270762750.172;
  double events_MC_tW_top = 277241050.839;
  double events_MC_W_jets = 29981320;
  double events_MC_dilep = 4984888995.13;

  double cross_MC_t_top = 136.02;
  double cross_MC_t_antitop = 80.95;
  double cross_MC_DY_50 = 6025.2;
  double cross_MC_DY_1050 = 22635.1;
  double cross_MC_semilep = 365.3;
  double cross_MC_tW_antitop = 35.9;
  double cross_MC_tW_top = 35.9;
  double cross_MC_W_jets = 61526.7;
  double cross_MC_dilep = 88.2;

  Histo_1->Scale(1.0/events_MC_t_top*cross_MC_t_top*lumi);
  Histo_2->Scale(1.0/events_MC_t_antitop*cross_MC_t_antitop*lumi);
  Histo_3->Scale(1.0/events_MC_DY_50*cross_MC_DY_50*lumi);
  Histo_4->Scale(1.0/events_MC_DY_1050*cross_MC_DY_1050*lumi);
  Histo_5->Scale(1.0/events_MC_semilep*cross_MC_semilep*lumi);
  Histo_6->Scale(1.0/events_MC_tW_antitop*cross_MC_tW_antitop*lumi);
  Histo_7->Scale(1.0/events_MC_tW_top*cross_MC_tW_top*lumi);
  Histo_8->Scale(1.0/events_MC_W_jets*cross_MC_W_jets*lumi);
  Histo_9->Scale(1.0/events_MC_dilep*cross_MC_dilep*lumi);


  TH1D *singletop = new TH1D ("stack single top", "", nbins, xmin, xmax );
  TH1D *ttbar_tw= new TH1D ("stack ttbar_tw", "", nbins, xmin, xmax );
  TH1D *QCD = new TH1D ("stack QCD", "", nbins, xmin, xmax );
  TH1D *w_z_gamma = new TH1D ("stack w_z_gamma", "", nbins, xmin, xmax );
  TH1D *all_others = new TH1D ("all_others", "", nbins, xmin, xmax );

  singletop->Add(Histo_1);
  singletop->Add(Histo_2);

  ttbar_tw->Add(Histo_5);
  ttbar_tw->Add(Histo_9);
  ttbar_tw->Add(Histo_6);
  ttbar_tw->Add(Histo_7);


  w_z_gamma->Add(Histo_3);
  w_z_gamma->Add(Histo_4);
  w_z_gamma->Add(Histo_8);

  QCD->Add(Histo_10);

  all_others->Add(singletop);
  all_others->Add(ttbar_tw);
  all_others->Add(w_z_gamma);
  double integ_all_others = all_others->Integral(1,nbins);
  double integ_QCD = QCD->Integral(1,nbins);
  double integ_data = DATA->Integral(1,nbins);

  double constrain1 = integ_all_others/integ_data;

  cout<<"integral QCD ="<<integ_QCD<<endl;

  cout<<"all others proportion = "<<integ_all_others/integ_data<<endl;
  cout<<"integral DATA ="<<integ_data<<endl;

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");


  double newscale0 = integ_data*param0/integ_all_others;
  double newscale1 = integ_data*param1/integ_QCD;
  cout<<"newscale0 = "<< newscale0<<endl;
  cout<<"newscale1 = "<< newscale1<<endl;
  singletop->Scale(newscale0);
  ttbar_tw->Scale(newscale0);
  w_z_gamma->Scale(newscale0);
  QCD->Scale(newscale1);

  double error_Wjets_tot[nbins];
  double error_ttbar_tw_tot[nbins];
  double error_tchan_tot[nbins];
  double error_multijet_tot[nbins];

  for (int i = 0 ; i< nbins ; i++)
  {

    error_Wjets_tot[i] = sqrt(pow(3200,2) + pow(2800,2))/(33400+30700)*w_z_gamma->GetBinContent(i+1);
    error_ttbar_tw_tot[i] = sqrt(pow(1400,2)+ pow(1500,2))/(84500+84800)*ttbar_tw->GetBinContent(i+1);
    error_tchan_tot[i] = sqrt(pow(820,2)+pow(880,2)+pow(3,2)+pow(2,2))/(17720+27+25+11460)*singletop->GetBinContent(i+1);
    error_multijet_tot[i] = 0.2*QCD->GetBinContent(i+1);

  }

  TH1D *all_error_up = new TH1D("error_up","error_up",nbins,xmin,xmax);
  TH1D *all_histo = new TH1D ("multijet", "multijet", nbins, xmin, xmax );
  TH1D *all_m1_histo = new TH1D ("w_z_gamma", "w_z_gamma", nbins, xmin, xmax );
  TH1D *all_m2_histo = new TH1D ("ttbar_tw", "ttbar_tw", nbins, xmin, xmax );
  TH1D *all_m3_histo = new TH1D ("singletop", "singletop", nbins, xmin, xmax );


  all_histo->Add(singletop);
  all_histo->Add(ttbar_tw);
  all_histo->Add(w_z_gamma);
  all_histo->Add(QCD);

  all_m1_histo->Add(singletop);
  all_m1_histo->Add(ttbar_tw);
  all_m1_histo->Add(w_z_gamma);

  all_m2_histo->Add(singletop);
  all_m2_histo->Add(ttbar_tw);

  all_m3_histo->Add(singletop);

  all_histo->SetLineColor(kBlack);
  all_histo->SetFillColor(kGray);

  all_m1_histo->SetLineColor(kBlack);
  all_m1_histo->SetFillColor(8);

  all_m2_histo->SetLineColor(kBlack);
  all_m2_histo->SetFillColor(kOrange);

  all_m3_histo->SetLineColor(kBlack);
  all_m3_histo->SetFillColor(kRed);

  all_error_up->SetFillStyle(3005);
  all_error_up->SetFillColor(kBlack);
  all_error_up->SetMarkerStyle(0);
  for (int i = 0 ; i<nbins ; i++)
  {
    all_error_up->SetBinContent(i+1, all_histo->GetBinContent(i+1));
    all_error_up->SetBinError(i+1,sqrt( pow(error_Wjets_tot[i],2) + pow(error_ttbar_tw_tot[i],2) + pow(error_tchan_tot[i],2) +pow(error_multijet_tot[i],2)));
  }

  DATA->SetLineWidth(2);
  DATA->SetLineColor(kBlack);
  all_histo->SetXTitle(legendX.c_str());
  all_histo->SetYTitle(legendY.c_str());
  double max = all_histo->GetMaximum();
  all_histo->SetAxisRange(0,max*1.3,"Y");
  if ( variable != "boosted_reco_top_mass")
  {
  all_histo->SetAxisRange(0,max*1.3,"Y");
  }

  all_histo->Draw("HIST");
  all_m1_histo->Draw("SAME HIST");
  all_m2_histo->Draw("SAME HIST");
  all_m3_histo->Draw("SAME HIST");
  all_error_up->Draw("E2 SAME");
  DATA->Draw("E SAME");

  string legendtitle1 = region+" Region";
  float lx0 = 0.65;
  float ly0 = 0.65;
  float lx1 = 0.99;
  float ly1 = 0.99;

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle1.c_str());
  string legendEntry1 = "DATA";
  string legendEntry2 = "Multijet (DATA)";
  string legendEntry3 = "W/Z/#gamma* + jets";
  string legendEntry4 = "t#bar{t}/tw";
  string legendEntry5 = "Single top t-channel";



  legend->SetFillColor(kWhite);
  legend->AddEntry(DATA->GetName(), legendEntry1.c_str(), "l");
  legend->AddEntry(all_histo->GetName(), legendEntry2.c_str(), "f");
  legend->AddEntry(all_m1_histo->GetName(), legendEntry3.c_str(), "f");
  legend->AddEntry(all_m2_histo->GetName(), legendEntry4.c_str(), "f");
  legend->AddEntry(all_m3_histo->GetName(), legendEntry5.c_str(), "f");
  legend->Draw("SAME");

  Canvas->Print(Name.c_str());


  string root_name = "./results/heppy/top_reco/ATCGRoot/"+variable+"_"+lepton+"_"+region+".root";

    TFile* file_histo = new TFile(root_name.c_str(),"RECREATE");

    file_histo->cd();
    DATA->SetStats(kTRUE);
    DATA->Write("data_obs");
    singletop->Write("sig");
    ttbar_tw->Write("ttbar_oth");
    w_z_gamma->Write("Wjets");
    QCD->Write("QCD");
    all_histo->Write("total");
    file_histo->Close();
    cout<<"integral = "<<DATA->Integral(1,nbins)<<endl;
    cout<<"nentry = "<<DATA->GetEntries()<<endl;
}

void Stacked_histo_reversed_Fit(TTree* t1, TTree* t2, TTree* t3, TTree* t4 , TTree* t5, TTree* t6, TTree* t7, TTree* t8, TTree* t9, TTree* t10, TTree* t11,
   string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string Name, string region){

  double lumi = 41500;
  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "Histo_2");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "Histo_3");
  TH1D* Histo_4 = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection, "Histo_4");
  TH1D* Histo_5 = GetHistoWeight(t5, variable, nbins, xmin, xmax, selection, "Histo_5");
  TH1D* Histo_6 = GetHistoWeight(t6, variable, nbins, xmin, xmax, selection, "Histo_6");
  TH1D* Histo_7 = GetHistoWeight(t7, variable, nbins, xmin, xmax, selection, "Histo_7");
  TH1D* Histo_8 = GetHistoWeight(t8, variable, nbins, xmin, xmax, selection, "Histo_8");
  TH1D* Histo_9 = GetHistoWeight(t9, variable, nbins, xmin, xmax, selection, "Histo_9");
  TH1D* Histo_10 = GetHistoWeight(t10, variable, nbins, xmin, xmax, selection, "Histo_10");
  TH1D* DATA = GetHistoWeight(t11, variable, nbins, xmin, xmax, selection, "Histo_29");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);
  Histo_4->SetStats(kFALSE);
  Histo_5->SetStats(kFALSE);
  Histo_6->SetStats(kFALSE);
  Histo_7->SetStats(kFALSE);
  Histo_8->SetStats(kFALSE);
  Histo_9->SetStats(kFALSE);
  Histo_10->SetStats(kFALSE);

  DATA->SetStats(kFALSE);

  double events_MC_t_top = 5982064;
  double events_MC_t_antitop = 3675910;
  double events_MC_DY_50 = 491250940213;
  double events_MC_DY_1050 = 315866;
  double events_MC_semilep = 32414370098;
  double events_MC_tW_antitop = 270762750.172;
  double events_MC_tW_top = 277241050.839;
  double events_MC_W_jets = 29981320;
  double events_MC_dilep = 4984888995.13;

  double cross_MC_t_top = 136.02;
  double cross_MC_t_antitop = 80.95;
  double cross_MC_DY_50 = 6025.2;
  double cross_MC_DY_1050 = 22635.1;
  double cross_MC_semilep = 365.3;
  double cross_MC_tW_antitop = 35.9;
  double cross_MC_tW_top = 35.9;
  double cross_MC_W_jets = 61526.7;
  double cross_MC_dilep = 88.2;

  Histo_1->Scale(1.0/events_MC_t_top*cross_MC_t_top*lumi);
  Histo_2->Scale(1.0/events_MC_t_antitop*cross_MC_t_antitop*lumi);
  Histo_3->Scale(1.0/events_MC_DY_50*cross_MC_DY_50*lumi);
  Histo_4->Scale(1.0/events_MC_DY_1050*cross_MC_DY_1050*lumi);
  Histo_5->Scale(1.0/events_MC_semilep*cross_MC_semilep*lumi);
  Histo_6->Scale(1.0/events_MC_tW_antitop*cross_MC_tW_antitop*lumi);
  Histo_7->Scale(1.0/events_MC_tW_top*cross_MC_tW_top*lumi);
  Histo_8->Scale(1/events_MC_W_jets*cross_MC_W_jets*lumi);
  Histo_9->Scale(1.0/events_MC_dilep*cross_MC_dilep*lumi);


  //TH1D *test = new TH1D ("stack single top", "", nbins, xmin, xmax );
  TH1D *singletop = new TH1D ("stack single top", "", nbins, xmin, xmax );
  TH1D *ttbar_tw= new TH1D ("stack ttbar_tw", "", nbins, xmin, xmax );
  TH1D *QCD = new TH1D ("stack QCD", "", nbins, xmin, xmax );
  TH1D *w_z_gamma = new TH1D ("stack w_z_gamma", "", nbins, xmin, xmax );
  TH1D *all_others = new TH1D ("all_others", "", nbins, xmin, xmax );

  singletop->Add(Histo_1);
  singletop->Add(Histo_2);

  ttbar_tw->Add(Histo_5);
  ttbar_tw->Add(Histo_9);
  ttbar_tw->Add(Histo_6);
  ttbar_tw->Add(Histo_7);


  w_z_gamma->Add(Histo_3);
  w_z_gamma->Add(Histo_4);
  w_z_gamma->Add(Histo_8);

  QCD->Add(Histo_10);
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  all_others->Add(singletop);
  all_others->Add(ttbar_tw);
  all_others->Add(w_z_gamma);

  double integ_all_others = all_others->Integral(1,nbins);
  double integ_QCD = QCD->Integral(1,nbins);
  double integ_data = DATA->Integral(1,nbins);

  double constrain1 = integ_all_others/integ_data;

  cout<<"integral QCD ="<<integ_QCD<<endl;

  cout<<"all others proportion = "<<integ_all_others/integ_data<<endl;
  cout<<"integral DATA ="<<integ_data<<endl;


  TObjArray *mc = new TObjArray(4);

  mc->Add(all_others);
  mc->Add(QCD);
  TFractionFitter* fit = new TFractionFitter(DATA,mc);


  fit->Constrain(0,constrain1-0.05,constrain1+0.05);
  fit->Constrain(1,0.0,0.5);


  fit->SetRangeX(1,nbins);
  Int_t status = fit->Fit();
  cout<<"fit status:"<<status<<endl;

  double param0, param1;
  double param0_error, param1_error;

  fit->GetResult(0,param0,param0_error);
  fit->GetResult(1,param1,param1_error);


  double newscale0 = integ_data*param0/integ_all_others;
  double newscale1 = integ_data*param1/integ_QCD;

  double full_error_0 = newscale0 * ( sqrt(integ_data)/integ_data + sqrt(integ_all_others)/integ_all_others + param0_error/param0 );
  double full_error_1 = newscale1 * ( sqrt(integ_data)/integ_data + sqrt(integ_QCD)/integ_QCD + param1_error/param1 );

  cout<<"newscale0="<<newscale0<<endl;
  cout<<"newscale0 error = "<<full_error_0<<endl;

  cout<<"newscale1="<<newscale1<<endl;
  cout<<"newscale1 error = "<<full_error_1<<endl;
  singletop->Scale(newscale0);
  ttbar_tw->Scale(newscale0);
  w_z_gamma->Scale(newscale0);
  QCD->Scale(newscale1);

  double error_Wjets_tot[nbins];
  double error_ttbar_tw_tot[nbins];
  double error_tchan_tot[nbins];
  double error_multijet_tot[nbins];

  for (int i = 0 ; i< nbins ; i++)
  {

    error_Wjets_tot[i] = sqrt(pow(3200,2) + pow(2800,2))/(33400+30700)*w_z_gamma->GetBinContent(i+1);
    error_ttbar_tw_tot[i] = sqrt(pow(1400,2)+ pow(1500,2))/(84500+84800)*ttbar_tw->GetBinContent(i+1);
    error_tchan_tot[i] = sqrt(pow(820,2)+pow(880,2)+pow(3,2)+pow(2,2))/(17720+27+25+11460)*singletop->GetBinContent(i+1);
    error_multijet_tot[i] = 0.2*QCD->GetBinContent(i+1);

  }
  TH1D *all_error_up = new TH1D("error_up","error_up",nbins,xmin,xmax);

  TH1D *all_histo = new TH1D ("multijet", "multijet", nbins, xmin, xmax );
  TH1D *all_m1_histo = new TH1D ("w_z_gamma", "w_z_gamma", nbins, xmin, xmax );
  TH1D *all_m2_histo = new TH1D ("ttbar_tw", "ttbar_tw", nbins, xmin, xmax );
  TH1D *all_m3_histo = new TH1D ("singletop", "singletop", nbins, xmin, xmax );


  all_histo->Add(singletop);
  all_histo->Add(ttbar_tw);
  all_histo->Add(w_z_gamma);
  all_histo->Add(QCD);

  all_m1_histo->Add(singletop);
  all_m1_histo->Add(ttbar_tw);
  all_m1_histo->Add(w_z_gamma);

  all_m2_histo->Add(singletop);
  all_m2_histo->Add(ttbar_tw);

  all_m3_histo->Add(singletop);

  all_histo->SetLineColor(kBlack);
  all_histo->SetFillColor(kGray);

  all_m1_histo->SetLineColor(kBlack);
  all_m1_histo->SetFillColor(8);

  all_m2_histo->SetLineColor(kBlack);
  all_m2_histo->SetFillColor(kOrange);

  all_m3_histo->SetLineColor(kBlack);
  all_m3_histo->SetFillColor(kRed);

  all_error_up->SetFillStyle(3005);
  all_error_up->SetFillColor(kBlack);
  all_error_up->SetMarkerStyle(0);

  for (int i = 0 ; i<nbins ; i++)
  {
    all_error_up->SetBinContent(i+1, all_histo->GetBinContent(i+1));
    all_error_up->SetBinError(i+1,sqrt( pow(error_Wjets_tot[i],2) + pow(error_ttbar_tw_tot[i],2) + pow(error_tchan_tot[i],2) +pow(error_multijet_tot[i],2)));
  }


  DATA->SetLineWidth(2);
  DATA->SetLineColor(kBlack);
  all_histo->SetXTitle(legendX.c_str());
  all_histo->SetYTitle(legendY.c_str());
  double max = all_histo->GetMaximum();
  all_histo->SetAxisRange(0,max*1.2,"Y");

  all_histo->Draw("HIST");
  all_m1_histo->Draw("SAME HIST");
  all_m2_histo->Draw("SAME HIST");
  all_m3_histo->Draw("SAME HIST");

  DATA->Draw("E SAME");
  all_error_up->Draw("E2 SAME");
  string legendtitle1 = region + " Region";
  float lx0 = 0.65;
  float ly0 = 0.65;
  float lx1 = 0.99;
  float ly1 = 0.99;

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle1.c_str());
  string legendEntry1 = "DATA";
  string legendEntry2 = "Multijet (DATA)";
  //string legendEntry3 = "all others MC";
  string legendEntry3 = "W/Z/#gamma* + jets";
  string legendEntry4 = "t#bar{t}/tw";
  string legendEntry5 = "Single top t-channel";



  legend->SetFillColor(kWhite);
  legend->AddEntry(DATA->GetName(), legendEntry1.c_str(), "l");
  legend->AddEntry(all_histo->GetName(), legendEntry2.c_str(), "f");
  legend->AddEntry(all_m1_histo->GetName(), legendEntry3.c_str(), "f");
  legend->AddEntry(all_m2_histo->GetName(), legendEntry4.c_str(), "f");
  legend->AddEntry(all_m3_histo->GetName(), legendEntry5.c_str(), "f");
  legend->Draw("SAME");

  Canvas->Print(Name.c_str());

}
void ATGCRoot(TTree* t1, TTree* t2, TTree* t3, TTree* t4 , TTree* t5, TTree* t6, TTree* t7, TTree* t8, TTree* t9, TTree* t10, TTree* t11,
   string variable, int nbins, double xmin, double xmax, string selection, double param0, double param1, string lepton, string name_file_variable , string EFT){

  double lumi = 41500;

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_4 = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_5 = GetHistoWeight(t5, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_6 = GetHistoWeight(t6, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_7 = GetHistoWeight(t7, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_8 = GetHistoWeight(t8, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_9 = GetHistoWeight(t9, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_10 = GetHistoWeight(t10, variable, nbins, xmin, xmax, selection, "");
  TH1D* DATA = GetHistoWeight(t11, variable, nbins, xmin, xmax, selection, "");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);
  Histo_4->SetStats(kFALSE);
  Histo_5->SetStats(kFALSE);
  Histo_6->SetStats(kFALSE);
  Histo_7->SetStats(kFALSE);
  Histo_8->SetStats(kFALSE);
  Histo_9->SetStats(kFALSE);
  Histo_10->SetStats(kFALSE);

  DATA->SetStats(kFALSE);

  double events_MC_t_top = 5982064;
  double events_MC_t_antitop = 3675910;
  double events_MC_DY_50 = 491250940213;
  double events_MC_DY_1050 = 315866;
  double events_MC_semilep = 32414370098;
  double events_MC_tW_antitop = 270762750.172;
  double events_MC_tW_top = 277241050.839;
  double events_MC_W_jets = 29981320;
  double events_MC_dilep = 4984888995.13;

  double cross_MC_t_top = 136.02;
  double cross_MC_t_antitop = 80.95;
  double cross_MC_DY_50 = 6025.2;
  double cross_MC_DY_1050 = 22635.1;
  double cross_MC_semilep = 365.3;
  double cross_MC_tW_antitop = 35.9;
  double cross_MC_tW_top = 35.9;
  double cross_MC_W_jets = 61526.7;
  double cross_MC_dilep = 88.2;

  Histo_1->Scale(1.0/events_MC_t_top*cross_MC_t_top*lumi);
  Histo_2->Scale(1.0/events_MC_t_antitop*cross_MC_t_antitop*lumi);
  Histo_3->Scale(1.0/events_MC_DY_50*cross_MC_DY_50*lumi);
  Histo_4->Scale(1.0/events_MC_DY_1050*cross_MC_DY_1050*lumi);
  Histo_5->Scale(1.0/events_MC_semilep*cross_MC_semilep*lumi);
  Histo_6->Scale(1.0/events_MC_tW_antitop*cross_MC_tW_antitop*lumi);
  Histo_7->Scale(1.0/events_MC_tW_top*cross_MC_tW_top*lumi);
  Histo_8->Scale(1.0/events_MC_W_jets*cross_MC_W_jets*lumi);
  Histo_9->Scale(1.0/events_MC_dilep*cross_MC_dilep*lumi);


  TH1D *singletop = new TH1D ("stack single top", "", nbins, xmin, xmax );
  TH1D *ttbar_tw= new TH1D ("stack ttbar_tw", "", nbins, xmin, xmax );
  TH1D *QCD = new TH1D ("stack QCD", "", nbins, xmin, xmax );
  TH1D *w_z_gamma = new TH1D ("stack w_z_gamma", "", nbins, xmin, xmax );
  TH1D *all_others = new TH1D ("all_others", "", nbins, xmin, xmax );
  TH1D *all_histo = new TH1D ("all_histo", "", nbins, xmin, xmax );
  TH1D *DATA_obs = new TH1D("data_obs","", nbins, xmin, xmax);

  singletop->Add(Histo_1);
  singletop->Add(Histo_2);

  ttbar_tw->Add(Histo_5);
  ttbar_tw->Add(Histo_9);
  ttbar_tw->Add(Histo_6);
  ttbar_tw->Add(Histo_7);


  w_z_gamma->Add(Histo_3);
  w_z_gamma->Add(Histo_4);
  w_z_gamma->Add(Histo_8);

  QCD->Add(Histo_10);
  DATA_obs->Add(DATA);

  all_others->Add(singletop);
  all_others->Add(ttbar_tw);
  all_others->Add(w_z_gamma);
  double integ_all_others = all_others->Integral(1,20);
  double integ_QCD = QCD->Integral(1,20);
  double integ_data = DATA->Integral(1,20);


  TCanvas* Canvas = new TCanvas("Canvas","Canvas");

  double newscale0 = integ_data*param0/integ_all_others;
  double newscale1 = integ_data*param1/integ_QCD;
  singletop->Scale(newscale0);
  ttbar_tw->Scale(newscale0);
  w_z_gamma->Scale(newscale0);
  QCD->Scale(newscale1);

  all_histo->Add(singletop);
  all_histo->Add(ttbar_tw);
  all_histo->Add(w_z_gamma);
  all_histo->Add(QCD);

  string root_name = "./results/heppy/top_reco/ATGCRoot/"+name_file_variable+"_"+EFT+".root";

    TFile* file_histo = new TFile(root_name.c_str(),"RECREATE");

    file_histo->cd();
    DATA->SetStats(kTRUE);
    DATA_obs->Write("data_obs");
    singletop->Write("diboson");
    ttbar_tw->Write("ttbar_oth");
    w_z_gamma->Write("Wjets");
    QCD->Write("QCD");
    all_histo->Write("total");
    file_histo->Close();
    cout<<"integral = "<<DATA->Integral(1,nbins)<<endl;
}

int main (){

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  //tdrStyle->SetEndErrorSize(2);
  //tdrStyle->SetErrorMarker(20);
  //tdrStyle->SetErrorX(0.);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.13);
  tdrStyle->SetPadRightMargin(0.05);

  // For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.04);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.05);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.03, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  tdrStyle->cd();

int dataset = 31;
string suffix[dataset];
suffix[0] = "MC_t_top";
suffix[1] = "MC_t_antitop";
suffix[2] = "MC_DY_50";
suffix[3] = "MC_DY_1050";
suffix[4] = "MC_semilep";
suffix[5] = "MC_tW_antitop";
suffix[6] = "MC_tW_top";
suffix[7] = "MC_W+Jets";
suffix[8] = "MC_dilep";
suffix[9] = "MC_QCD_EM_20to30";
suffix[10] = "MC_QCD_EM_30to50";
suffix[11] = "MC_QCD_EM_50to80";
suffix[12] = "MC_QCD_EM_80to120";
suffix[13] = "MC_QCD_EM_120to170";
suffix[14] = "MC_QCD_EM_170to300";
suffix[15] = "MC_QCD_EM_300toInf";
suffix[16] = "MC_QCD_Mu_15to20";
suffix[17] = "MC_QCD_Mu_20to30";
suffix[18] = "MC_QCD_Mu_30to50";
suffix[19] = "MC_QCD_Mu_50to80";
suffix[20] = "MC_QCD_Mu_80to120";
suffix[21] = "MC_QCD_Mu_120to170";
suffix[22] = "MC_QCD_Mu_170to300";
suffix[23] = "MC_QCD_Mu_300to470";
suffix[24] = "MC_QCD_Mu_470to600";
suffix[25] = "MC_QCD_Mu_600to800";
suffix[26] = "MC_QCD_Mu_800to1000";
suffix[27] = "MC_QCD_Mu_1000toInf";
suffix[28] = "DATA_region";
suffix[29] = "multijet_elec_0_muon_iso";
suffix[30] = "multijet_muon_0_elec_iso";
/*
suffix[0] = "MC_t_top_QCD";
suffix[1] = "MC_t_antitop_QCD";
suffix[2] = "MC_DY_50_QCD";
suffix[3] = "MC_DY_1050_QCD";
suffix[4] = "MC_semilep_QCD";
suffix[5] = "MC_tW_antitop_QCD";
suffix[6] = "MC_tW_top_QCD";
suffix[7] = "MC_W+Jets_QCD";
suffix[8] = "MC_dilep_QCD";
suffix[9] = "MC_QCD_EM_20to30_QCD";
suffix[10] = "MC_QCD_EM_30to50_QCD";
suffix[11] = "MC_QCD_EM_50to80_QCD";
suffix[12] = "MC_QCD_EM_80to120_QCD";
suffix[13] = "MC_QCD_EM_120to170_QCD";
suffix[14] = "MC_QCD_EM_170to300_QCD";
suffix[15] = "MC_QCD_EM_300toInf_QCD";
suffix[16] = "MC_QCD_Mu_15to20_QCD";
suffix[17] = "MC_QCD_Mu_20to30_QCD";
suffix[18] = "MC_QCD_Mu_30to50_QCD";
suffix[19] = "MC_QCD_Mu_50to80_QCD";
suffix[20] = "MC_QCD_Mu_80to120_QCD";
suffix[21] = "MC_QCD_Mu_120to170_QCD";
suffix[22] = "MC_QCD_Mu_170to300_QCD";
suffix[23] = "MC_QCD_Mu_300to470_QCD";
suffix[24] = "MC_QCD_Mu_470to600_QCD";
suffix[25] = "MC_QCD_Mu_600to800_QCD";
suffix[26] = "MC_QCD_Mu_800to1000_QCD";
suffix[27] = "MC_QCD_Mu_1000toInf_QCD";
suffix[28] = "DATA_QCD";*/

	TFile* fInput[dataset];
	TTree* tInput[dataset];
	string inputName;

  for (int i=0; i<dataset; i++)
    {
      cout<<"treating "<<suffix[i]<<endl;
      inputName = "data/heppy/output/output_" + suffix[i] +".root";
      fInput[i] = new TFile(inputName.c_str(),"READ");
      tInput[i] = (TTree*) fInput[i]->Get("events");
    }

//--------------W + jets control region-------------------------//

/*  Compare_1Histos(tInput[0], "region_0_jet1_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","jet 1 pt in w+jets control region","results/heppy/"+suffix[0]+"_w+jets_j1_pt.pdf");
  Compare_1Histos(tInput[0], "region_0_jet2_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","jet 2 pt in w+jets control region","results/heppy/"+suffix[0]+"_w+jets_j2_pt.pdf");
  Compare_1Histos(tInput[0], "region_0_jet1_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","jet 1 eta in w+jets control region","results/heppy/"+suffix[0]+"_w+jets_j1_eta.pdf");
  Compare_1Histos(tInput[0], "region_0_jet2_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","jet 2 eta in w+jets control region","results/heppy/"+suffix[0]+"_w+jets_j2_eta.pdf");
  Compare_1Histos(tInput[0], "region_0_elec_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","elec pt in w+jets control region","results/heppy/"+suffix[0]+"_w+jets_elec_pt.pdf");
  Compare_1Histos(tInput[0], "region_0_elec_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","elec eta in w+jets control region","results/heppy/"+suffix[0]+"_w+jets_elec_eta.pdf");
  Compare_1Histos(tInput[0], "region_0_muon_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","muon pt in w+jets control region","results/heppy/"+suffix[0]+"_w+jets_muon_pt.pdf");
  Compare_1Histos(tInput[0], "region_0_muon_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","muon eta in w+jets control region","results/heppy/"+suffix[0]+"_w+jets_muon_eta.pdf");*/
  /*Compare_1Histos(tInput[0], "region_0_mtw", 100, 0, 400, "weight_one_lepton", "mtw", "normalized", "legendUpRight", "MC","W transverse mass in w+jets control region","results/heppy/"+suffix[0]+"_w+jets_mtw.pdf");
  Compare_1Histos(tInput[0], "region_0_elec_iso", 100, -0.001, 0.06, "weight_one_lepton", "electron isolation", "normalized", "legendUpRight", "MC","electron relative isolation in w+jets control region","results/heppy/"+suffix[0]+"_w+jets_elec_iso.pdf");
  Compare_1Histos(tInput[0], "region_0_muon_iso", 100, -0.001, 0.15, "weight_one_lepton", "muon isolation", "normalized", "legendUpRight", "MC","muon relative isolation in w+jets control region","results/heppy/"+suffix[0]+"_w+jets_muon_iso.pdf");*/

//-------------Signal region-------------------------------------//

  /*Compare_1Histos(tInput[0], "region_1_jet1_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","jet 1 pt in signal region","results/heppy/"+suffix[0]+"_signal_j1_pt.pdf");
  Compare_1Histos(tInput[0], "region_1_jet2_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","jet 2 pt in signal region","results/heppy/"+suffix[0]+"_signal_j2_pt.pdf");
  Compare_1Histos(tInput[0], "region_1_b1_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC", "b-tagged jet 1 pt in signal region","results/heppy/"+suffix[0]+"_signal_b1_pt.pdf");

  Compare_1Histos(tInput[0], "region_1_jet1_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","jet 1 eta in signal region","results/heppy/"+suffix[0]+"_signal_j1_eta.pdf");
  Compare_1Histos(tInput[0], "region_1_jet2_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","jet 2 eta in signal region","results/heppy/"+suffix[0]+"_signal_j2_eta.pdf");
  Compare_1Histos(tInput[0], "region_1_b1_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC", "b-tagged jet 1 eta in signal region","results/heppy/"+suffix[0]+"_signal_b1_eta.pdf");

  Compare_1Histos(tInput[0], "region_1_elec_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","elec pt in signal region","results/heppy/"+suffix[0]+"_signal_elec_pt.pdf");
  Compare_1Histos(tInput[0], "region_1_elec_eta", 20, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","elec eta in signal region","results/heppy/"+suffix[0]+"_signal_elec_eta.pdf");

  Compare_1Histos(tInput[0], "region_1_muon_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","muon pt in signal region","results/heppy/"+suffix[0]+"_signal_muon_pt.pdf");
  Compare_1Histos(tInput[0], "region_1_muon_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","muon eta in signal region","results/heppy/"+suffix[0]+"_signal_muon_eta.pdf");*/
  /*Compare_1Histos(tInput[0], "region_1_mtw", 100, 0, 400, "weight_one_lepton", "mtw", "normalized", "legendUpRight", "MC","W transverse mass in signal region","results/heppy/"+suffix[0]+"_signal_mtw.pdf");
  Compare_1Histos(tInput[0], "region_1_elec_iso", 100, -0.001, 0.06, "weight_one_lepton", "electron isolation", "normalized", "legendUpRight", "MC","electron relative isolation in signal region","results/heppy/"+suffix[0]+"_signal_elec_iso.pdf");
  Compare_1Histos(tInput[0], "region_1_muon_iso", 100, -0.001, 0.15, "weight_one_lepton", "muon isolation", "normalized", "legendUpRight", "MC","muon relative isolation in signal region","results/heppy/"+suffix[0]+"_signal_muon_iso.pdf");*/
  /*Compare_1Histos(tInput[0], "final_region", 3, 0,3, "final_natureLepton==1", "region for elec", "normalized", "legendUpRight", "DATA/MC","region","results/heppy/"+suffix[0]+"_elec.pdf");
  Compare_1Histos(tInput[0], "final_region", 3, 0,3, "final_natureLepton==2", "region for muon", "normalized", "legendUpRight", "DATA/MC","region","results/heppy/"+suffix[0]+"_muon.pdf");
  Compare_1Histos(tInput[0], "final_region", 3, 0,3, "1", "region all leptons", "normalized", "legendUpRight", "DATA/MC","region","results/heppy/"+suffix[0]+"_leptons.pdf");*/

//---------------ttbar control region---------------------------//

  /*Compare_1Histos(tInput[0], "region_2_jet1_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","jet 1 pt in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_j1_pt.pdf");
  Compare_1Histos(tInput[0], "region_2_jet2_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","jet 2 pt in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_j2_pt.pdf");
  Compare_1Histos(tInput[0], "region_2_jet3_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","jet 3 pt in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_j3_pt.pdf");

  Compare_1Histos(tInput[0], "region_2_jet1_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","jet 1 eta in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_j1_eta.pdf");
  Compare_1Histos(tInput[0], "region_2_jet2_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","jet 2 eta in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_j2_eta.pdf");
  Compare_1Histos(tInput[0], "region_2_jet3_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","jet 3 eta in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_j3_eta.pdf");

  Compare_1Histos(tInput[0], "region_2_b1_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","b-tagged jet 1 pt in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_b1_pt.pdf");
  Compare_1Histos(tInput[0], "region_2_b2_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","b-tagged jet 2 pt in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_b2_pt.pdf");

  Compare_1Histos(tInput[0], "region_2_b1_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","b-tagged jet 1 eta in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_b1_eta.pdf");
  Compare_1Histos(tInput[0], "region_2_b2_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","b-tagged jet 2 eta in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_b2_eta.pdf");

  Compare_1Histos(tInput[0], "region_2_elec_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","elec pt in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_elec_pt.pdf");
  Compare_1Histos(tInput[0], "region_2_elec_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","elec eta in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_elec_eta.pdf");

  Compare_1Histos(tInput[0], "region_2_muon_pt", 100, 0, 400, "weight_one_lepton", "pt", "normalized", "legendUpRight", "MC","muon pt in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_muon_pt.pdf");
  Compare_1Histos(tInput[0], "region_2_muon_eta", 40, 0, 4, "weight_one_lepton", "eta", "normalized", "legendUpRight", "MC","muon eta in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_muon_eta.pdf");*/
  /*Compare_1Histos(tInput[0], "region_2_mtw", 100, 0, 400, "weight_one_lepton", "mtw", "normalized", "legendUpRight", "MC","W transverse mass in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_mtw.pdf");
  Compare_1Histos(tInput[0], "region_2_elec_iso", 100, -0.001, 0.06, "weight_one_lepton", "electron isolation", "normalized", "legendUpRight", "MC","electron relative isolation in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_elec_iso.pdf");
  Compare_1Histos(tInput[0], "region_2_muon_iso", 100, -0.001, 0.15, "weight_one_lepton", "muon isolation", "normalized", "legendUpRight", "MC","muon relative isolation in ttbar control region","results/heppy/"+suffix[0]+"_ttbar_muon_iso.pdf");*/
  //---------------Data------------------//

  //--------------W + jets control region-------------------------//

  /*Compare_1Histos(tInput[1], "region_0_jet1_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA","jet 1 pt in w+jets control region","results/heppy/"+suffix[1]+"_w+jets_j1_pt.pdf");
  Compare_1Histos(tInput[1], "region_0_jet2_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA","jet 2 pt in w+jets control region","results/heppy/"+suffix[1]+"_w+jets_j2_pt.pdf");
  Compare_1Histos(tInput[1], "region_0_jet1_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","jet 1 eta in w+jets control region","results/heppy/"+suffix[1]+"_w+jets_j1_eta.pdf");
  Compare_1Histos(tInput[1], "region_0_jet2_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","jet 2 eta in w+jets control region","results/heppy/"+suffix[1]+"_w+jets_j2_eta.pdf");
  Compare_1Histos(tInput[1], "region_0_elec_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA","elec pt in w+jets control region","results/heppy/"+suffix[1]+"_w+jets_elec_pt.pdf");
  Compare_1Histos(tInput[1], "region_0_elec_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","elec eta in w+jets control region","results/heppy/"+suffix[1]+"_w+jets_elec_eta.pdf");
  Compare_1Histos(tInput[1], "region_0_muon_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA","muon pt in w+jets control region","results/heppy/"+suffix[1]+"_w+jets_muon_pt.pdf");
  Compare_1Histos(tInput[1], "region_0_muon_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","muon eta in w+jets control region","results/heppy/"+suffix[1]+"_w+jets_muon_eta.pdf");*/
  /*Compare_1Histos(tInput[1], "region_0_mtw", 100, 0, 400, "1", "mtw", "normalized", "legendUpRight", "DATA","W transverse mass in w+jets control region","results/heppy/"+suffix[1]+"_w+jets_mtw.pdf");
  Compare_1Histos(tInput[1], "region_0_elec_iso", 100, -0.001, 0.06, "1", "electron isolation", "normalized", "legendUpRight", "DATA","electron relative isolation in w+jets control region","results/heppy/"+suffix[1]+"_w+jets_elec_iso.pdf");
  Compare_1Histos(tInput[1], "region_0_muon_iso", 100, -0.001, 0.15, "1", "muon isolation", "normalized", "legendUpRight", "DATA","muon relative isolation in w+jets control region","results/heppy/"+suffix[1]+"_w+jets_muon_iso.pdf");*/
//-------------Signal region-------------------------------------//

  /*Compare_1Histos(tInput[1], "region_1_jet1_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA","jet 1 pt in signal region","results/heppy/"+suffix[1]+"_signal_j1_pt.pdf");
  Compare_1Histos(tInput[1], "region_1_jet2_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA","jet 2 pt in signal region","results/heppy/"+suffix[1]+"_signal_j2_pt.pdf");
  Compare_1Histos(tInput[1], "region_1_b1_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA", "b-tagged jet 1 pt in signal region","results/heppy/"+suffix[1]+"_signal_b1_pt.pdf");

  Compare_1Histos(tInput[1], "region_1_jet1_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","jet 1 eta in signal region","results/heppy/"+suffix[1]+"_signal_j1_eta.pdf");
  Compare_1Histos(tInput[1], "region_1_jet2_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","jet 2 eta in signal region","results/heppy/"+suffix[1]+"_signal_j2_eta.pdf");
  Compare_1Histos(tInput[1], "region_1_b1_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA", "b-tagged jet 1 eta in signal region","results/heppy/"+suffix[1]+"_signal_b1_eta.pdf");

  Compare_1Histos(tInput[1], "region_1_elec_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA","elec pt in signal region","results/heppy/"+suffix[1]+"_signal_elec_pt.pdf");
  Compare_1Histos(tInput[1], "region_1_elec_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","elec eta in signal region","results/heppy/"+suffix[1]+"_signal_elec_eta.pdf");

  Compare_1Histos(tInput[1], "region_1_muon_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA","muon pt in signal region","results/heppy/"+suffix[1]+"_signal_muon_pt.pdf");
  Compare_1Histos(tInput[1], "region_1_muon_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","muon eta in signal region","results/heppy/"+suffix[1]+"_signal_muon_eta.pdf");*/
  /*Compare_1Histos(tInput[1], "region_1_mtw", 100, 0, 400, "1", "mtw", "normalized", "legendUpRight", "DATA","W transverse mass in signal region","results/heppy/"+suffix[1]+"_signal_mtw.pdf");
  Compare_1Histos(tInput[1], "region_1_elec_iso", 100, -0.001, 0.06, "1", "electron isolation", "normalized", "legendUpRight", "DATA","electron relative isolation in signal region","results/heppy/"+suffix[1]+"_signal_elec_iso.pdf");
  Compare_1Histos(tInput[1], "region_1_muon_iso", 100, -0.001, 0.15, "1", "muon isolation", "normalized", "legendUpRight", "DATA","muon relative isolation in signal region","results/heppy/"+suffix[1]+"_signal_muon_iso.pdf");*/

//---------------------ttbar control region---------------------------//

  /*Compare_1Histos(tInput[1], "region_2_jet1_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA","jet 1 pt in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_j1_pt.pdf");
  Compare_1Histos(tInput[1], "region_2_jet2_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA","jet 2 pt in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_j2_pt.pdf");
  Compare_1Histos(tInput[1], "region_2_jet3_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA","jet 3 pt in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_j3_pt.pdf");

  Compare_1Histos(tInput[1], "region_2_jet1_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","jet 1 eta in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_j1_eta.pdf");
  Compare_1Histos(tInput[1], "region_2_jet2_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","jet 2 eta in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_j2_eta.pdf");
  Compare_1Histos(tInput[1], "region_2_jet3_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","jet 3 eta in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_j3_eta.pdf");

  Compare_1Histos(tInput[1], "region_2_b1_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA","b-tagged jet 1 pt in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_b1_pt.pdf");
  Compare_1Histos(tInput[1], "region_2_b2_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA","b-tagged jet 2 pt in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_b2_pt.pdf");

  Compare_1Histos(tInput[1], "region_2_b1_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","b-tagged jet 1 eta in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_b1_eta.pdf");
  Compare_1Histos(tInput[1], "region_2_b2_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","b-tagged jet 2 eta in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_b2_eta.pdf");

  Compare_1Histos(tInput[1], "region_2_elec_pt", 100, 0, 40, "1", "pt", "normalized", "legendUpRight", "DATA","elec pt in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_elec_pt.pdf");
  Compare_1Histos(tInput[1], "region_2_elec_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","elec eta in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_elec_eta.pdf");

  Compare_1Histos(tInput[1], "region_2_muon_pt", 100, 0, 400, "1", "pt", "normalized", "legendUpRight", "DATA","muon pt in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_muon_pt.pdf");
  Compare_1Histos(tInput[1], "region_2_muon_eta", 40, 0, 4, "1", "eta", "normalized", "legendUpRight", "DATA","muon eta in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_muon_eta.pdf");*/
  /*Compare_1Histos(tInput[1], "region_2_mtw", 100, 0, 400, "1", "mtw", "normalized", "legendUpRight", "DATA","W transverse mass in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_mtw.pdf");
  Compare_1Histos(tInput[1], "region_2_elec_iso", 100, -0.001, 0.06, "1", "electron isolation", "normalized", "legendUpRight", "DATA","electron relative isolation in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_elec_iso.pdf");
  Compare_1Histos(tInput[1], "region_2_muon_iso", 100, -0.001, 0.15, "1", "muon isolation", "normalized", "legendUpRight", "DATA","muon relative isolation in ttbar control region","results/heppy/"+suffix[1]+"_ttbar_muon_iso.pdf");*/

  //--------------------DATA driven multijet background----------------------//

  //------------------------------Muon--------------------------------------//

  /*Compare_1Histos(tInput[2], "region_0_mtw", 100, 0, 400, "1", "mtw", "normalized", "legendUpRight", "DATA with reversed cut on muon Irel > 0.2","W transverse mass in w+jets control region","results/heppy/"+suffix[2]+"_w+jets_mtw.pdf");
  Compare_1Histos(tInput[2], "region_0_muon_iso", 500, 0.19, 70, "1", "muon isolation", "normalized", "legendUpRight", "DATA","muon relative isolation in w+jets control region","results/heppy/"+suffix[2]+"_w+jets_muon_iso.pdf");

  //-----------------------------Electron---------------------------------//

  Compare_1Histos(tInput[3], "region_0_mtw", 100, 0, 400, "1", "mtw", "normalized", "legendUpRight", "DATA with events failing the veto for electron ID and Irel < 0.85","W transverse mass in w+jets control region","results/heppy/"+suffix[3]+"_w+jets_mtw.pdf");
  Compare_1Histos(tInput[3], "region_0_elec_iso", 100, -0.001, 0.85, "1", "electron isolation", "normalized", "legendUpRight", "DATA","electron relative isolation in w+jets control region","results/heppy/"+suffix[3]+"_w+jets_elec_iso.pdf");*/


  //-----------------------------DATA+MC---------------------------------//

/*Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
  , tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_region", 3, 0, 3, "final_natureLepton==1", "region elec", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_region_elec.pdf");
*/
/*Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_mtw", 20, 0, 300, "final_natureLepton == 2 && final_region==0", "M_{t,W} electron", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_mtw_muon.pdf");
Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_elec_pt", 15, 0, 200, "final_natureLepton==1 && final_region==1", "p_{T} electron (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_pt_elec.pdf",1,1,"signal");
/*Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_elec_eta", 20, 0, 4, "final_region==1 && final_natureLepton==1", "|#eta| electron", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_eta_elec.pdf");

Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_jet1_eta", 20, 0, 4, "final_region==1 && final_natureLepton ==1", "|#eta| jet electron", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_eta_jet1_elec.pdf");
Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_jet1_pt", 20, 0, 300, "final_region==1 && final_natureLepton ==1", "p_T jet electron", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_pt_jet1_elec.pdf");

Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_b1_pt", 20, 0, 300, "final_region==1 && final_natureLepton ==1", "p_T jet electron", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_pt_bjet1_elec.pdf");

Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_b1_eta", 20, 0, 4, "final_region==1 && final_natureLepton ==1", "|#eta| b-jet electron", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_eta_bjet1_elec.pdf");


Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
  , tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_region", 3, 0, 3, "final_natureLepton==2", "region muon", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_region_muon.pdf");

Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_mtw", 20, 0, 300, "final_region==1 && final_natureLepton == 2 ", "M_{t,W} muon", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_mtw_muon.pdf");
Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_muon_pt", 20, 0, 300, "final_region==1 && final_natureLepton==2 ", "p_{T} muon (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_pt_muon.pdf",1,1,"signal");
/*Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_muon_eta", 20, 0, 4, "final_region==1 && final_natureLepton==2", "|#eta| muon", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_eta_muon.pdf");
*//*
Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_jet1_eta", 20, 0, 4, "final_region==1 && final_natureLepton ==2", "|#eta| jet muon", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_eta_jet1_muon.pdf");
Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_jet1_pt", 20, 0, 300, "final_region==1 && final_natureLepton ==2", "p_T jet muon", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_pt_jet1_muon.pdf");

Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_b1_pt", 20, 0, 300, "final_region==1 && final_natureLepton ==2", "p_T b-jet muon", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_pt_bjet1_muon.pdf");

Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_b1_eta", 20, 0, 4, "final_region==1 && final_natureLepton ==2", "|#eta| b-jet muon", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_eta_bjet1_muon.pdf");*/
/*
Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_jet1_eta", 20, -5, 5, "final_region==0 && final_natureLepton ==2", "|#eta|", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_eta_jet1_muon_region0.pdf");
Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_jet2_eta", 20, -5, 5, "final_region==0 && final_natureLepton ==2", "|#eta|", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_eta_jet2_muon_region0.pdf");


//------------------------------------Comparison for DATA driven multijet background-------------------------//
/*Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "M_T_W", 20, 0, 300, "final_natureLepton == 1 && final_region==0", "M_{t,W} (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_mtw_elec_region0_DATA_reversed_postFit.pdf");

//----------------------------DATA + MC with multijet from DATA with Chi square fit----------------------///

/*Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "final_mtw", 20, 0, 300, "final_natureLepton == 1 && final_region==0", "M_{t,W} (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_mtw_elec_with_not_region0.pdf");*/

/*Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "M_T_W", 20, 0, 300, "final_natureLepton == 1 && final_region==0 && M_T_W < 600", "M_{t,W} (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_mtw_elec_without_not_region0_postFit.pdf");*/

/*Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "boosted_reco_top_mass", 20, 0, 500, "final_natureLepton == 1 && final_region==1 && 0 < boosted_reco_top_mass < 600", "Masse top (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_top_mass_elec_with_not_region1_TFrac.pdf");
*/

//---------------------------------variable with electrons--------------------//
/*Stacked_histo_reversed_Fit(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "M_T_W", 20, 0, 300, "final_natureLepton == 1 && final_region==0 && M_T_W < 600", "M_{W,T} (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_mtw_elec_error_postFit.pdf", "W+Jets");

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "M_T_W", 20, 0, 300, "final_natureLepton == 1 && final_region==0 ", "M_{W,T} (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_mtw_elec_region0_postFit.pdf", 0.800042, 0.202856, "elec", "W+Jets");

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "jet_not_b_eta", 20, -4.5, 4.5, "final_natureLepton == 1 && final_region==1 ", "#eta", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_not_b_eta_elec_region1_postFit.pdf", 0.800042, 0.202856, "elec", "Signal");
Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "M_T_top", 20, 0, 400, "final_natureLepton == 1 && final_region==2 ", "M_{top,T} (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_mt_top_elec_region2_postFit.pdf", 0.800042, 0.202856, "elec", "ttbar");
/*Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "M_T_tot", 20, 0, 150, "final_natureLepton == 1 && final_region==2 ", "M_{tot,T} (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_mt_tot_elec_region2_postFit.pdf", 0.800042, 0.201963, "elec", "ttbar");


Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "boosted_reco_cosTheta", 20, -1, 1, "final_natureLepton == 1 && final_region==1 && boosted_reco_cosTheta <=1 ", "cos(#theta)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_cosTheta_elec_without_not_region1_postFit.pdf", 0.800042, 0.201963, "elec", "Signal");

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "boosted_reco_PhiStar", 20, 0, 6.2831, "final_natureLepton == 1 && final_region==1", "#phi^{*}", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_PhiStar_elec_without_not_region1_postFit.pdf", 0.800042, 0.201963, "elec" , "Signal");

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "boosted_reco_cosThetaStar", 20, -1, 1, "final_natureLepton == 1 && final_region==1 && boosted_reco_cosThetaStar <=1 ", "cos(#theta^{*})", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_cosThetaStar_elec_without_not_region1_postFit.pdf",  0.800042, 0.201963, "elec", "Signal");

/*Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "boosted_reco_top_mass", 20, 0, 500, "final_natureLepton == 1 && final_region==1 && 0 < boosted_reco_top_mass < 600", "Masse top (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_top_mass_elec_without_not_region1_postFit.pdf",  0.800042, 0.201963, "elec", "signal");

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "final_elec_pt", 15, 0, 200, "final_natureLepton==1 && final_region==1", "p_{T} electron (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_pt_elec_postfit.pdf",  0.800042, 0.202856, "elec", "Signal");
/*
Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "Delta_R", 10, 0, 6, "final_natureLepton==1 && final_region==1", "#DeltaR", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_Delta_R_elec_region1_postfit.pdf",  0.800042, 0.201963, "elec", "Signal ");

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "final_jet1_eta", 10, -5, 5, "final_natureLepton==1 && final_region==0", "|#eta|", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_jet1_eta_elec_region0_postfit.pdf",  0.800042, 0.201963, "elec", "W+Jets");
Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "final_jet2_eta", 10, -5, 5, "final_natureLepton==1 && final_region==0", "|#eta|", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_jet2_eta_elec_region0_postfit.pdf",  0.800042, 0.201963, "elec", "W+Jets");

/*

ATGCRoot(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "boosted_reco_cosTheta", 20, -1, 1, "final_natureLepton == 1 && final_region==1 && boosted_reco_cosTheta <=1 ", 0.800042, 0.201963, "elec", "cosTheta", "ctwi");

ATGCRoot(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "boosted_reco_PhiStar", 5,  0, 6.2831, "final_natureLepton == 1 && final_region==1 ", 0.800042, 0.201963, "elec", "PhiStar", "ctwi");

//---------------------------------variable with muons--------------------//*/
/*Stacked_histo_reversed_Fit(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "M_T_W", 20, 0, 300, "final_natureLepton == 2 && final_region==0 && M_T_W < 600", "M_{W,T} (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_mtw_muon_error_postFit.pdf", "W+Jets");

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "M_T_W", 20, 0, 300, "final_natureLepton == 2 && final_region==0 ", "M_{W,T} (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_mtw_muon_region0_postFit.pdf", 0.741288,0.263463, "muon", "W+Jets");
Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "jet_not_b_eta", 20, -5, 5, "final_natureLepton == 2 && final_region==1 ", "#eta", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_not_b_eta_muon_region1_postFit.pdf", 0.741288,0.263463, "muon", "Signal");
Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "M_T_top", 20, 0, 400, "final_natureLepton == 2 && final_region==2 ", "M_{top,T} (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_mt_top_muon_region2_postFit.pdf", 0.741288,0.263463, "muon", "ttbar");
Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "M_T_tot", 20, 0, 140, "final_natureLepton == 2 && final_region==2 ", "M_{tot,T} (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_mt_tot_muon_region2_postFit.pdf", 0.741288,0.263463, "muon", "ttbar");
/*Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "Delta_R", 10, 0, 6, "final_natureLepton == 2 && final_region==1 ", "#DeltaR", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_Delta_R_muon_region1_postFit.pdf", 0.741288,0.263463, "muon", "Signal region");
Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "final_jet1_pt", 50, 0, 300, "final_region == 2 && final_natureLepton == 2 && abs(final_jet1_eta) <3.0 && abs(final_jet1_eta)>2.7 ", "p_{T}", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_jet1_pt_muon_region2_postFit.pdf", 0.741288,0.263463, "muon", "ttbar");
Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "final_jet2_pt", 50, 0, 300, "final_region == 2 && final_natureLepton == 2 && abs(final_jet2_eta) <3.0 && abs(final_jet2_eta)>2.7 ", "p_{T}", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_jet2_pt_muon_region2_postFit.pdf", 0.741288,0.263463, "muon", "ttbar");
Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "final_jet3_pt", 50, 0, 300, "final_region == 2 && final_natureLepton == 2 && abs(final_jet3_eta) <3.0 && abs(final_jet3_eta)>2.7 ", "p_{T}", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_jet3_pt_muon_region2_postFit.pdf", 0.741288,0.263463, "muon", "ttbar");

/*Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "final_jet2_eta", 10, -5, 5, "final_region == 0 && final_natureLepton == 2 ", "|#eta|", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_jet2_eta_muon_region0_postFit.pdf", 0.741288,0.263463, "muon", "WJets");
Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "final_jet3_eta", 10, -5, 5, "final_region == 1 && final_natureLepton == 2 ", "|#eta|", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_jet3_eta_muon_region1_postFit.pdf", 0.741288,0.263463, "muon", "signal");

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "boosted_reco_cosTheta", 20, -1, 1, "final_natureLepton == 2 && final_region==1 && boosted_reco_cosTheta <=1 ", "cos(#theta)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_cosTheta_muon_with_not_region1_postFit.pdf", 0.741288,0.263463,"muon", "Signal");


Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "boosted_reco_PhiStar", 20, 0, 6.2831, "final_natureLepton == 2 && final_region==1", "#phi^{*}", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_PhiStar_muon_with_not_region1_postFit.pdf", 0.741288,0.263463, "muon", "Signal");

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "boosted_reco_cosThetaStar", 20, -1, 1, "final_natureLepton == 2 && final_region==1 && boosted_reco_cosThetaStar <=1 ", "cos(#theta^{*})", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_cosThetaStar_muon_with_not_region1_postFit.pdf", 0.741288,0.263463, "muon", "Signal");

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "boosted_reco_top_mass", 20, 0, 500, "final_natureLepton == 2 && final_region==1 && 0 < boosted_reco_top_mass < 600", "Masse top (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_top_mass_muon_with_not_region1_postFit.pdf", 0.741288,0.263463, "muon","signal");
/*Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "final_natureLepton", 1, 2, 2, "final_natureLepton == 2 && final_region==1 && charge_lepton == -1", "nature lepton", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_number_muon_with_not_region1_postFit.pdf", 0.651319,0.349910);

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "final_natureLepton", 1, 2, 2, "final_natureLepton == 2 && final_region==1 && charge_lepton == 1", "nature lepton", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_number_antimuon_with_not_region1_postFit.pdf", 0.651319,0.349910);
*/
    //---------------------with MC--------------------------//


/*Stacked_histo_Fit(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
    , tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "M_T_W", 20, 0, 300, "final_natureLepton == 1 && final_region==0", "M_{t,W} (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_multijetMC_mtw_elec_region0_MC_postFit.pdf");
/*
Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "boosted_reco_cosTheta", 20, -1, 1, "final_region==1 && final_natureLepton == 1 && boosted_reco_cosTheta < 2", "cos(#theta)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_multijetMC_cosTheta_elec_region1_postFit.pdf",0.7516,0.236);

Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "boosted_reco_top_mass", 20, 0, 500, "final_region==1 && final_natureLepton == 1 && 0 < boosted_reco_top_mass < 500 ", "M_{top}", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_multijetMC_top_mass_elec_region1_postFit.pdf",0.7516,0.236);
/*
Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "boosted_reco_cosThetaStar", 20, -1, 1, "final_region==1 && final_natureLepton == 1 && boosted_reco_cosThetaStar < 2", "cos(#theta)^*", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_multijetMC_cosThetaStar_elec_region1_postFit.pdf",0.7516,0.236);
*//*
Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_jet1_eta", 20, -5, 5, " final_region == 2 && final_natureLepton == 2 ", "|#eta|", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_multijetMC_jet1_eta_muon_region2_postFit.pdf",0.860288,0.139786,"ttbar");
Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_jet2_eta", 20, -5, 5, " final_region == 2 && final_natureLepton == 2 ", "|#eta|", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_multijetMC_jet2_eta_muon_region2_postFit.pdf",0.860288,0.139786, "ttbar");
Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "final_jet3_eta", 20, -5, 5, " final_region == 2 && final_natureLepton == 2 ", "|#eta|", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_multijetMC_jet3_eta_muon_region2_postFit.pdf",0.860288,0.139786, "ttbar");
/*
//------------------------------Variable with muons------------------------------//
Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "boosted_reco_cosTheta", 20, -1, 1, "final_natureLepton == 2 && final_region==1 && -2< boosted_reco_cosTheta < 2 && final_iso_muon > 0.2 ", "cos(#theta)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_cosTheta_muon_with_not_region1_TFrac.pdf");

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "boosted_reco_top_mass", 20, 0, 500, "final_natureLepton == 2 && final_region==1 && 0 < boosted_reco_top_mass < 600 final_iso_muon > 0.2", "Masse top (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_top_mass_muon_with_not_region1_TFrac.pdf");

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[30],
tInput[28], "reco_top_pt", 20, 0, 300, "final_natureLepton == 2 && final_region==1 && reco_top_pt <300 final_iso_muon > 0.2 ", "Masse top (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_top_pt_muon_with_not_region1_TFrac.pdf");



Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "boosted_reco_cosTheta", 20, -1, 1, "final_natureLepton == 2 && final_region==1 && boosted_reco_cosTheta <=1 ", "cos(#theta)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_cosTheta_muon_with_not_region1.pdf");

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "boosted_reco_top_mass", 20, 0, 500, "final_natureLepton == 2 && final_region==1 && 0 < boosted_reco_top_mass < 600 ", "Masse top (GeV)", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_top_mass_muon_with_not_region1.pdf");

Stacked_histo_reversed(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29],
tInput[28], "reco_top_pt", 20, 0, 500, "final_natureLepton == 2 && final_region==1 && 0 < reco_top_pt < 600", "Masse top (GeV) final_iso_muon > 0.2", "number of events", "legendUpRight", "DATA/MC","results/heppy/top_reco/DATA_MC_top_pt_muon_with_not_region1.pdf");
*/

//-------------------------------FULL QCD--------------------------------------//

/*Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15], tInput[16], tInput[17]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "log(final_muon_iso)", 50, -10, 10, "final_region==0 && final_natureLepton==2", "log I_{#mu,rel}", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_QCD_Iso_muon.pdf",1,1,"W+Jets");
/*Stacked_histo(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[9], tInput[10], tInput[11], tInput[12], tInput[13], tInput[14], tInput[15], tInput[16], tInput[17]
, tInput[18], tInput[19], tInput[20], tInput[21], tInput[22], tInput[23], tInput[24], tInput[25], tInput[26], tInput[27], tInput[28], "log(final_elec_iso)", 50, -10, 10, "final_region==0 && final_natureLepton==1", "log I_{elec,rel}", "number of events", "legendUpRight", "DATA/MC","results/heppy/DATA_MC_QCD_Iso_elec.pdf",1,1,"W+Jets");
*/


//--------------------------- Chi Square--------------------------------------//

int number_of_bins = 20;
string variable = "boosted_reco_PhiStar";
string name_file_plot = "signal_proc_PhiStar_ctwi";
double xminimum =0;
double xmaximum= 6.2831;
/*string variable = "boosted_reco_cosThetaStar";
string name_file_plot = "signal_proc_cosThetaStar_cbwi";
double xminimum =-1;
double xmaximum= 1;*/
/*
double syst_lumi = 0.023;
double syst_singletop = 0.06;
double syst_ttbar_tw = 0.02;
double syst_w_z_gamma = 0.08;
double syst_QCD = 0.10;

*/
double stat = normalizeWithChi(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29], tInput[30],
tInput[28],variable,number_of_bins,xminimum,xmaximum,"final_region==1",0.800042, 0.201963, 0.741288,0.262616, 1.0, 1.0, 1.0, 1.0, 1.0, name_file_plot, "cos(#theta^{*})","number of events","results/heppy/DATA_MC_signal_EFT_PhiStar.pdf", "signal" );
/*
double statsyst1[5];
double statsyst2[5];

statsyst1[0] = normalizeWithChi(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29], tInput[30],
tInput[28],variable,number_of_bins,xminimum,xmaximum,"final_region==1",0.800042, 0.201963, 0.741288,0.262616,1.0+syst_lumi, 1.0, 1.0, 1.0, 1.0, name_file_plot, "ratio","#Chi^{2}");
statsyst1[1] = normalizeWithChi(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29], tInput[30],
tInput[28],variable,number_of_bins,xminimum,xmaximum,"final_region==1 ",0.800042, 0.201963, 0.741288,0.262616 ,1.0, 1.0+syst_singletop, 1.0, 1.0, 1.0, name_file_plot, "ratio","#Chi^{2}");
statsyst1[2] = normalizeWithChi(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29], tInput[30],
tInput[28],variable,number_of_bins,xminimum,xmaximum,"final_region==1 ",0.800042, 0.201963, 0.741288,0.262616 ,1.0, 1.0, 1.0+syst_ttbar_tw, 1.0, 1.0, name_file_plot, "ratio","#Chi^{2}");
statsyst1[3] = normalizeWithChi(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29], tInput[30],
tInput[28],variable,number_of_bins,xminimum,xmaximum,"final_region==1 ",0.800042, 0.201963, 0.741288,0.262616 ,1.0, 1.0, 1.0, 1.0+syst_w_z_gamma, 1.0, name_file_plot, "ratio","#Chi^{2}");
statsyst1[4] = normalizeWithChi(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29], tInput[30],
tInput[28],variable,number_of_bins,xminimum,xmaximum,"final_region==1 ",0.800042, 0.201963, 0.741288,0.262616,1.0, 1.0, 1.0, 1.0, 1.0+syst_QCD, name_file_plot, "ratio","#Chi^{2}");

statsyst2[0] = normalizeWithChi(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29], tInput[30],
tInput[28],variable,number_of_bins,xminimum,xmaximum,"final_region==1 ",0.800042, 0.201963, 0.741288,0.262616 ,1.0-syst_lumi, 1.0, 1.0, 1.0, 1.0, name_file_plot, "ratio","#Chi^{2}");
statsyst2[1] = normalizeWithChi(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29], tInput[30],
tInput[28],variable,number_of_bins,xminimum,xmaximum,"final_region==1 ",0.800042, 0.201963, 0.741288,0.262616 ,1.0, 1.0-syst_singletop, 1.0, 1.0, 1.0, name_file_plot, "ratio","#Chi^{2}");
statsyst2[2] = normalizeWithChi(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29], tInput[30],
tInput[28],variable,number_of_bins,xminimum,xmaximum,"final_region==1 ",0.800042, 0.201963, 0.741288,0.262616, 1.0, 1.0, 1.0-syst_ttbar_tw, 1.0, 1.0, name_file_plot, "ratio","#Chi^{2}");
statsyst2[3] = normalizeWithChi(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29], tInput[30],
tInput[28],variable,number_of_bins,xminimum,xmaximum,"final_region==1 ",0.800042, 0.201963, 0.741288,0.262616 ,1.0, 1.0, 1.0, 1.0-syst_w_z_gamma, 1.0, name_file_plot, "ratio","#Chi^{2}");
statsyst2[4] = normalizeWithChi(tInput[0], tInput[1], tInput[2], tInput[3], tInput[4], tInput[5], tInput[6], tInput[7], tInput[8], tInput[29], tInput[30],
tInput[28],variable,number_of_bins,xminimum,xmaximum,"final_region==1 ",0.800042, 0.201963, 0.741288,0.262616 ,1.0, 1.0, 1.0, 1.0,  1.0-syst_QCD, name_file_plot, "ratio","#Chi^{2}");

double syst[5];
double max_statsyst[5];
double statsyst_tot=stat*stat;
double squaredSum=0;

for (int i = 0; i <5;i++ )
{
   max_statsyst[i] = sqrt(abs(statsyst1[i]*statsyst1[i])+abs(statsyst2[i]*statsyst2[i]))/2;
   syst[i] = sqrt(abs(max_statsyst[i]*max_statsyst[i] - stat*stat));
   statsyst_tot += syst[i]*syst[i];
   squaredSum += syst[i]*syst[i];
}

cout<<"stat mu ="<<stat<<endl;
cout<<"syst ="<<sqrt(squaredSum)<<endl;
cout<<"tot ="<<sqrt(statsyst_tot)<<endl;*/

  return 0;

}
