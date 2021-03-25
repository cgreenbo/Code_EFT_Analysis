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
#include <string>
#include <TGraph.h>
#include <TF1.h>
#include <fstream>
#include <TSystem.h>
using namespace std;

bool EFT_vs_SM = false;
bool reweight = true;

bool normalization = true; //Turn on normalization to 1
bool log_fity = false; //Turn on log scale in Y axis
bool diff_xsection = false; //Change this to have histograms normalized to Xsection
bool flows = true; //Change this to Active/Deactivate Overflow and Underflow for all Histograms


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

        string cutnew = "1 * (" + cut + ")";
        //string cutnew = "(" + cut + ")";

//

        t->Draw(variablenew.c_str(), cutnew.c_str());
        TH1D *histo = (TH1D*)gDirectory->Get("h");
  if(flows)
  {
		if (histo->GetEntries()==0) return histo;
  
		double underflow = histo->GetBinContent(0);
		//cout << "underflow="<<underflow<<endl;
		double val = 0;
		if (underflow>0) {
			val = histo->GetBinContent(1);
      //cout<<"val= "<<val<<endl;
			histo->SetBinContent(1, val+underflow);
			 histo->SetBinContent(0, 0);
		}
		double overflow = histo->GetBinContent(nbins+1);
    //cout<<"overflow= "<<overflow<<endl;
		if (overflow>0) {
		  val = histo->GetBinContent(nbins);
		  histo->SetBinContent(nbins+1, 0);
		  histo->SetBinContent(nbins, val+overflow);
		}
  }
  //cout << "Area="<<histo->Integral()<<endl;
	//cout << "Nevents="<<histo->GetEntries()<<endl;
        histo->SetName(name.c_str());
        histo->SetTitle(name.c_str());

        return histo;
}

void Ratio_EFT_SM(TTree* t1, TTree* t2, TTree* t3, TTree* t4, TTree* t5, string variable, string EFT, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string file_name, string lepton)
{
  TFile* file_output = new TFile(file_name.c_str(),"RECREATE");
  TTree* tree_file = new TTree("events","events");

  tree_file->Branch("");

  TH1D* Histo_SM = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_p1 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_p2 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_m1= GetHistoWeight(t4, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_m2 = GetHistoWeight(t5, variable, nbins, xmin, xmax, selection, "");

  TH1D* Ratio_p2_SM = (TH1D*)Histo_EFT_p2->Clone("Ratio_p2_SM");
  TH1D* Ratio_p1_SM = (TH1D*)Histo_EFT_p1->Clone("Ratio_p1_SM");
  TH1D* Ratio_m1_SM = (TH1D*)Histo_EFT_m1->Clone("Ratio_m1_SM");
  TH1D* Ratio_m2_SM = (TH1D*)Histo_EFT_m2->Clone("Ratio_m2_SM");

  Ratio_p2_SM->Divide(Histo_SM);
  Ratio_p1_SM->Divide(Histo_SM);
  Ratio_m2_SM->Divide(Histo_SM);
  Ratio_m1_SM->Divide(Histo_SM);

  bool div_Xsection = true;
  if(div_Xsection)
  {
    //int N=35.2863;
    int N = 1;
    Ratio_p2_SM->Scale(37.3893/N);
    Ratio_p1_SM->Scale(35.8126/N);
    Ratio_m2_SM->Scale(37.3838/N);
    Ratio_m1_SM->Scale(35.8029/N);
  }

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");

  Ratio_p2_SM->SetLineColor(kRed);
  Ratio_p2_SM->SetLineWidth(2);

  Ratio_p1_SM->SetLineColor(kBlue);
  Ratio_p1_SM->SetLineWidth(2);

  Ratio_m2_SM->SetLineColor(kOrange);
  Ratio_m2_SM->SetLineWidth(2);

  Ratio_m1_SM->SetLineColor(kBlack);
  Ratio_m1_SM->SetLineWidth(2);


  Ratio_p2_SM->SetXTitle(legendX.c_str());
  Ratio_p2_SM->SetYTitle(legendY.c_str());

  double max = (Ratio_p2_SM->GetMaximum()>Ratio_p1_SM->GetMaximum()) ? Ratio_p2_SM->GetMaximum() : Ratio_p1_SM->GetMaximum();

  Ratio_p2_SM->SetAxisRange((2-max)*1.1, max*1.25, "Y");
  Ratio_p2_SM->Draw("");
  Ratio_p1_SM->Draw("SAME");
  Ratio_m1_SM->Draw("SAME");
  Ratio_m2_SM->Draw("SAME");

  double lx0 = 0.6;
  double ly0 = 0.6;
  double lx1 = 0.99;
  double ly1 = 0.99;
  string legendtitle = "Value of the EFT";

  string eft_p2_legend = EFT + "/#\Lambda^{2} = 2 (TeV^{-2})";
  string eft_p1_legend = EFT + "/#\Lambda^{2} = 1 (TeV^{-2})";
  string eft_m1_legend = EFT + "/#\Lambda^{2} = -1 (TeV^{-2})";
  string eft_m2_legend = EFT + "/#\Lambda^{2} = -2 (TeV^{-2})";


  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.05);
  legend->AddEntry(Ratio_p2_SM->GetName(), eft_p2_legend.c_str(), "l");
  legend->AddEntry(Ratio_p1_SM->GetName(), eft_p1_legend.c_str(), "l");
  legend->AddEntry(Ratio_m1_SM->GetName(), eft_m1_legend.c_str(), "l");
  legend->AddEntry(Ratio_m2_SM->GetName(), eft_m2_legend.c_str(), "l");

  legend->Draw("SAME");

  TGraph** ratio_Histo = new TGraph*[nbins];

    //string TFile_name = "/results/ratio_madgraph/signal_proc_"+variable;
  TFile* ratio_file = new TFile(("signal_proc_"+variable+"_"+EFT+"_"+lepton+"_TF1.root").c_str(),"RECREATE");
  ratio_file->cd();
  Canvas->Print(file_name.c_str());

  int bin1 = 0;
  int bin2 = 1.25;


  for (int i = 1 ; i<=nbins ; i++)
  {
    //if(i == 3) continue;
    string number_plot = to_string(i); //Conversion of an int into a stream
    string file_name_eft = file_name;
    file_name_eft.insert(file_name.size()-4,"_"+number_plot);
    Canvas->Clear();
    //Canvas->SetLogy();
    ratio_Histo[i] = new TGraph(5);
    TF1* ratio_formula = new TF1(("bin_content_par1_"+number_plot).c_str(),"[0]+[1]*x+[2]*x*x", -3 , 3);
    
    bool Get_bins = false;
    if(Get_bins)
    {
      ratio_Histo[i]->SetPoint(0,-2,Histo_EFT_m2->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(1,-1,Histo_EFT_m1->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(2,0,Histo_SM->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(3,1,Histo_EFT_p1->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(4,2,Histo_EFT_p2->GetBinContent(i)/Histo_SM->GetBinContent(i));
    }
    else
    { Histo_SM->Divide(Histo_SM);
      Histo_SM->Scale(35.2863);
      ratio_Histo[i]->SetPoint(0,-2,Ratio_m2_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(1,-1,Ratio_m1_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(2,0,Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(3,1,Ratio_p1_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(4,2,Ratio_p2_SM->GetBinContent(i));
    }
    ratio_Histo[i]->SetMarkerStyle(kStar);

    tree_file->Fill();
    ratio_Histo[i]->Fit(ratio_formula);
    ratio_Histo[i]->Draw();
    //Canvas->Print(("test.pdf"+number_plot).c_str());
    
    /*
    if(i==1 || i==2 || i==3)
    { 
      TLegend* legend2 = new TLegend(0.6, 0.7, 0.89, 0.89, "");
      legend2->SetTextSize(0.05);
    
      if(i == 1) legend2->AddEntry(ratio_Histo[i]->GetName(), "0 < #phi^{*} < 1.25" ,"l");
      if(i == 2) legend2->AddEntry(ratio_Histo[i]->GetName(), "1.25 < #phi^{*} < 2.5" ,"l");
      if(i == 3) legend2->AddEntry(ratio_Histo[i]->GetName(), "2.5 < #phi^{*} < 4" ,"l");
      legend2->Draw(""); 
    }
    
    
    if(i==4 || i==5)
    {
      TLegend* legend2 = new TLegend(0.1, 0.2, 0.11, 0.11, "");
      if(i == 4) legend2->AddEntry(ratio_Histo[i]->GetName(), "4 < #phi^{*} < 5.5" ,"r");
      if(i == 5) legend2->AddEntry(ratio_Histo[i]->GetName(), "5.5 < #phi^{*} < 6.28" ,"r");
      legend2->Draw("");    
    }
    */

      TLegend* legend2 = new TLegend(0.6, 0.7, 0.89, 0.89, "");
      legend2->SetTextSize(0.05);
    
      //if(i == 1) legend2->AddEntry(ratio_Histo[i]->GetName(), "0 < #phi^{*} < 1.25" ,"r");
      //if(i == 2) legend2->AddEntry(ratio_Histo[i]->GetName(), "1.25 < #phi^{*} < 2.5" ,"r");
      //if(i == 3) legend2->AddEntry(ratio_Histo[i]->GetName(), "2.5 < #phi^{*} < 4" ,"r");
      //if(i == 4) legend2->AddEntry(ratio_Histo[i]->GetName(), "4 < #phi^{*} < 5.5" ,"r");
      //if(i == 5) legend2->AddEntry(ratio_Histo[i]->GetName(), "5.5 < #phi^{*} < 6.28" ,"r");
      //legend2->Draw("SAME");    

    
    ratio_Histo[i]->GetYaxis()->SetTitle("#sigma [pb]");
    ratio_Histo[i]->GetXaxis()->SetTitle("cbwi");



    
     

    
    Canvas->Print(file_name_eft.c_str());

    ratio_Histo[i]->Write(("bin_content_par1_"+number_plot).c_str());
    //ratio_formula->Write();
  }


  //file_output->Write();
  ratio_file->Close();




}

void Compare_3Histos(TTree* t1, TTree* t2, TTree* t3, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string legendEntry3, string legendEntry4, string Name){
  
  Name += variable;

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "Histo_2");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "Histo_3");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);

  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");

  if(log_fity) 
  {
    Canvas->SetLogy();
    Name += "_log"; 
    legendY += " log";
  }
  if(normalization) 
  {
    Double_t integral1 = 1/Histo_1->Integral();
    Histo_1->Scale(integral1);
    
    Double_t integral2 = 1/Histo_2->Integral();
    Histo_2->Scale(integral2);
    
    Double_t integral3 = 1/Histo_3->Integral();
    Histo_3->Scale(integral3);

    Name += "_normalized";
    legendY += " normalized";
  }

  if(diff_xsection)
  {
    //int N1 = Histo_1->GetEntries() , N2 = Histo_2->GetEntries() , N3 = Histo_3->GetEntries();
    int N_initial = 1000000;
    double scale1 = 35.2863/N_initial , scale2 = 37.3838/N_initial , scale3 = 37.3893/N_initial ;
    Histo_1->Scale(scale1);
    Histo_2->Scale(scale2);
    Histo_3->Scale(scale3);

    Name += "_xsection" ;
    legendY += " #sigma";
  }

  if (legendX == "Top Mass (GeV)")
    {
      if (normalization == false) Histo_1->SetAxisRange(0.0001,max*1.5,"Y");
    }
  if (legendX == "W Mass (GeV)")
  {
    Histo_1->SetAxisRange(0.001,max*1.5,"Y");
  }

  if(log_fity==false)
  {
    if(variable == "PhiStar")
    {
      Histo_1->SetMinimum(0);
      double max_phi = Histo_1->GetMaximum();
      Histo_1->SetMaximum(max_phi*1.5);
    }
    if(variable == "cosThetaStar")
    {
      Histo_1->SetMinimum(0);
      double max_cos = Histo_1->GetMaximum();
      Histo_1->SetMaximum(max_cos*1.28);
    }
  }

  else
  {
    if(variable == "PhiStar")
    {
      Histo_1->SetMinimum(0.1);
      double max_1 = Histo_1->GetMaximum();
      Histo_1->SetMaximum(max_1*40);

    }
    if(variable =="cosThetaStar")
    {
      Histo_1->SetMinimum(0.1);     
      double max_1 = Histo_1->GetMaximum();
      Histo_1->SetMaximum(max_1*40);
    }
  }
  


  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  
  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

  
  Histo_3->SetLineColor(kOrange);
  Histo_3->SetLineWidth(2);
  Histo_3->Draw("SAME");


 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.7;
	 lx1 = 0.6;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.6;
	 ly0 = 0.75; //Lower Y
	 lx1 = 0.99;
	 ly1 = 0.99; //Upper Y
  }
  
  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.04);
  if(legendEntry4 == "cbWi")
  {
    legend->AddEntry(Histo_1->GetName(), "C_{bW}^{I}/#\Lambda^{2} = 0 (TeV^{-2})","l");
    legend->AddEntry(Histo_2->GetName(), "C_{bW}^{I}/#\Lambda^{2} = -2 (TeV^{-2})","l");
    legend->AddEntry(Histo_3->GetName(), "C_{bW}^{I}/#\Lambda^{2} = 2 (TeV^{-2})","l");
  } 

  if(legendEntry4 == "ctWi")
  {
    legend->AddEntry(Histo_1->GetName(), "C_{tW}^{I}/#\Lambda^{2} = 0 (TeV^{-2})","l");
    legend->AddEntry(Histo_2->GetName(), "C_{tW}^{I}/#\Lambda^{2} = -2 (TeV^{-2})","l");
    legend->AddEntry(Histo_3->GetName(), "C_{tW}^{I}/#\Lambda^{2} = 2 (TeV^{-2})","l");
  } 
  //legend->SetLegendSize(0.5);

  legend->Draw("SAME");

  Name += ".gif";
  Canvas->Print(Name.c_str());
}

void Compare_1Histos(TTree* t1, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  Histo_1->SetStats(kFALSE);

  double max = Histo_1->GetMaximum();
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
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
  legend->AddEntry(Histo_1->GetName(), "MadGraph EFT = 0", "l");
  legend->AddEntry(Histo_2->GetName(), "MadSpin EFT = 0", "l");
  legend->AddEntry(Histo_3->GetName(), "MadGraph EFT = 2", "l");
  legend->AddEntry(Histo_4->GetName(), "MadSpin EFT = 2", "l");
  legend->Draw("SAME");

  Canvas->Print(Name.c_str());

  cout << "Histo1 mean: "<<Histo_1->GetMean()<<endl;
  cout << "Histo2 mean: "<<Histo_2->GetMean()<<endl;
  if (c>0) cout << "Histo3 mean: "<<Histo_3->GetMean()<<endl;

}

//New Fcts
inline bool exists_test0 (const std::string& name) 
{
  ifstream f(name.c_str());
  return f.good();
}


void Ratio_EFT_SM_7pts(TTree* t1, TTree* t2, TTree* t3, TTree* t4, TTree* t5, TTree* t6, TTree* t7, string variable, string EFT, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string file_name, string lepton)
{
  //Fcts to Plot EFT/SM for cbwi/ctwi = [-5;5]

  TFile* file_output = new TFile(file_name.c_str(),"RECREATE");
  TTree* tree_file = new TTree("events","events");

  tree_file->Branch("");

  TH1D* Histo_SM = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_p1 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_p2 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_m1= GetHistoWeight(t4, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_m2 = GetHistoWeight(t5, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_p5 = GetHistoWeight(t6, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_m5 = GetHistoWeight(t7, variable, nbins, xmin, xmax, selection, "");

  TH1D* Ratio_p2_SM = (TH1D*)Histo_EFT_p2->Clone("Ratio_p2_SM");
  TH1D* Ratio_p1_SM = (TH1D*)Histo_EFT_p1->Clone("Ratio_p1_SM");
  TH1D* Ratio_m1_SM = (TH1D*)Histo_EFT_m1->Clone("Ratio_m1_SM");
  TH1D* Ratio_m2_SM = (TH1D*)Histo_EFT_m2->Clone("Ratio_m2_SM");
  TH1D* Ratio_p5_SM = (TH1D*)Histo_EFT_p5->Clone("Ratio_p5_SM");
  TH1D* Ratio_m5_SM = (TH1D*)Histo_EFT_m5->Clone("Ratio_m5_SM");


  //Ratio_p2_SM->Divide(Histo_SM);
  //Ratio_p1_SM->Divide(Histo_SM);
  //Ratio_m2_SM->Divide(Histo_SM);
  //Ratio_m1_SM->Divide(Histo_SM);
  //Ratio_p5_SM->Divide(Histo_SM);
  //Ratio_m5_SM->Divide(Histo_SM);

  bool div_Xsection = true;
  if(div_Xsection)
  {
    double Xs_n5, Xs_n2, Xs_n1 ,Xs_SM = 35.2863;
    double Xs_p1, Xs_p2, Xs_p5;

    if(EFT == "cbwi")
    {
      Xs_n5 = 48.5315;
      Xs_n2 = 37.3699;
      Xs_n1 = 35.8251;
      Xs_p1 = 35.8091;
      Xs_p2 = 37.3955;
      Xs_p5 = 48.5184;
    }
    if(EFT == "ctwi")
    {
      Xs_n5 = 48.5683;
      Xs_n2 = 37.3838;
      Xs_n1 = 35.8029;
      Xs_p1 = 35.8126;
      Xs_p2 = 37.3893;
      Xs_p5 = 48.5945;
    }
    int N = 1000000;
    Ratio_p2_SM->Scale(Xs_p2/N);
    Ratio_p1_SM->Scale(Xs_p1/N);
    Ratio_m2_SM->Scale(Xs_n2/N);
    Ratio_m1_SM->Scale(Xs_n1/N);
    Ratio_p5_SM->Scale(Xs_p5/N);
    Ratio_m5_SM->Scale(Xs_n5/N);
    Histo_SM->Scale(Xs_SM/N);
  }

  bool ratio = true;
  if(ratio)
  {
    Ratio_p2_SM->Divide(Histo_SM);
    Ratio_p1_SM->Divide(Histo_SM);
    Ratio_m2_SM->Divide(Histo_SM);
    Ratio_m1_SM->Divide(Histo_SM);
    Ratio_p5_SM->Divide(Histo_SM);
    Ratio_m5_SM->Divide(Histo_SM);
    Histo_SM->Divide(Histo_SM);
  }
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");

  Ratio_p2_SM->SetLineColor(kRed);
  Ratio_p2_SM->SetLineWidth(2);

  Ratio_p1_SM->SetLineColor(kBlue);
  Ratio_p1_SM->SetLineWidth(2);

  Ratio_m2_SM->SetLineColor(kOrange);
  Ratio_m2_SM->SetLineWidth(2);

  Ratio_m1_SM->SetLineColor(kBlack);
  Ratio_m1_SM->SetLineWidth(2);

  Ratio_p5_SM->SetLineColor(kGreen);
  Ratio_p5_SM->SetLineWidth(2);

  Ratio_m5_SM->SetLineColor(kViolet);
  Ratio_m5_SM->SetLineWidth(2);


  Ratio_p2_SM->SetXTitle(legendX.c_str());
  Ratio_p2_SM->SetYTitle(legendY.c_str());

  double max = (Ratio_p5_SM->GetMaximum()>Ratio_m5_SM->GetMaximum()) ? Ratio_p5_SM->GetMaximum() : Ratio_m5_SM->GetMaximum();

  Ratio_p2_SM->SetAxisRange((2-max)*1.1, max*1.25, "Y");
  Ratio_p2_SM->Draw("");
  Ratio_p1_SM->Draw("SAME");
  Ratio_m1_SM->Draw("SAME");
  Ratio_m2_SM->Draw("SAME");
  Ratio_p5_SM->Draw("SAME");
  Ratio_m5_SM->Draw("SAME");

  double lx0 = 0.6;
  double ly0 = 0.6;
  double lx1 = 0.99;
  double ly1 = 0.99;
  string legendtitle = "Value of the EFT";

  string eft_p2_legend = EFT + "/#\Lambda^{2} = 2 (TeV^{-2})";
  string eft_p1_legend = EFT + "/#\Lambda^{2} = 1 (TeV^{-2})";
  string eft_m1_legend = EFT + "/#\Lambda^{2} = -1 (TeV^{-2})";
  string eft_m2_legend = EFT + "/#\Lambda^{2} = -2 (TeV^{-2})";
  string eft_p5_legend = EFT + "/#\Lambda^{2} = 5 (TeV^{-2})";
  string eft_m5_legend = EFT + "/#\Lambda^{2} = -5 (TeV^{-2})";



  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.05);
  legend->AddEntry(Ratio_p2_SM->GetName(), eft_p2_legend.c_str(), "l");
  legend->AddEntry(Ratio_p1_SM->GetName(), eft_p1_legend.c_str(), "l");
  legend->AddEntry(Ratio_m1_SM->GetName(), eft_m1_legend.c_str(), "l");
  legend->AddEntry(Ratio_m2_SM->GetName(), eft_m2_legend.c_str(), "l");
  legend->AddEntry(Ratio_p5_SM->GetName(), eft_p5_legend.c_str(), "l");
  legend->AddEntry(Ratio_m5_SM->GetName(), eft_m5_legend.c_str(), "l");

  legend->Draw("SAME");

  TGraph** ratio_Histo = new TGraph*[nbins];

    //string TFile_name = "/results/ratio_madgraph/signal_proc_"+variable;
  TFile* ratio_file = new TFile(("signal_proc_"+variable+"_"+EFT+"_"+lepton+"_TF1.root").c_str(),"RECREATE");
  ratio_file->cd();
  Canvas->Print(file_name.c_str());


  for (int i = 1 ; i<=nbins ; i++)
  {
    string number_plot = to_string(i); //Conversion of an int into a stream
    string file_name_eft = file_name;
    file_name_eft.insert(file_name.size()-4,"_"+number_plot);
    Canvas->Clear();
    //Canvas->SetLogy();
    ratio_Histo[i] = new TGraph(5);
    TF1* ratio_formula = new TF1(("bin_content_par1_"+number_plot).c_str(),"[0]+[1]*x+[2]*x*x", -3 , 3);
    
    bool Get_bins = false;
    if(Get_bins)
    {
      ratio_Histo[i]->SetPoint(0,-2,Histo_EFT_m2->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(1,-1,Histo_EFT_m1->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(2,0,Histo_SM->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(3,1,Histo_EFT_p1->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(4,2,Histo_EFT_p2->GetBinContent(i)/Histo_SM->GetBinContent(i));
    }
    else
    { 
      ratio_Histo[i]->SetPoint(0,-5,Ratio_m5_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(1,-2,Ratio_m2_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(2,-1,Ratio_m1_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(3,0,Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(4,1,Ratio_p1_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(5,2,Ratio_p2_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(6,5,Ratio_p5_SM->GetBinContent(i));
    }
    ratio_Histo[i]->SetMarkerStyle(kStar);

    tree_file->Fill();
    ratio_Histo[i]->Fit(ratio_formula);
    ratio_Histo[i]->Draw();
    //Canvas->Print(("test.pdf"+number_plot).c_str());

      TLegend* legend2 = new TLegend(0.49, 0.73, 0.75, 0.86, "");
      legend2->SetTextSize(0.035);
    if(variable == "PhiStar")
    {
      if(i == 1) legend2->AddEntry(ratio_Histo[i]->GetName(), "0 < #phi^{*} < 1.25" ,"r");
      if(i == 2) legend2->AddEntry(ratio_Histo[i]->GetName(), "1.25 < #phi^{*} < 2.5" ,"r");
      if(i == 3) legend2->AddEntry(ratio_Histo[i]->GetName(), "2.5 < #phi^{*} < 4" ,"r");
      if(i == 4) legend2->AddEntry(ratio_Histo[i]->GetName(), "4 < #phi^{*} < 5.5" ,"r");
      if(i == 5) legend2->AddEntry(ratio_Histo[i]->GetName(), "5.5 < #phi^{*} < 6.28" ,"r");
    }
    
    if(variable == "cosThetaStar")
    {
      if(i == 1) legend2->AddEntry(ratio_Histo[i]->GetName(), "-1 < cos(#theta^{*}) < -0.6" ,"r");
      if(i == 2) legend2->AddEntry(ratio_Histo[i]->GetName(), "-0.6 < cos(#theta^{*}) < -0.2" ,"r");
      if(i == 3) legend2->AddEntry(ratio_Histo[i]->GetName(), "-0.2 < cos(#theta^{*}) < 0.2" ,"r");
      if(i == 4) legend2->AddEntry(ratio_Histo[i]->GetName(), "0.2 < cos(#theta^{*}) < 0.6" ,"r");
      if(i == 5) legend2->AddEntry(ratio_Histo[i]->GetName(), "0.6 < cos(#theta^{*}) < 1" ,"r");

    }
    
    legend2->Draw("SAME");    

    ratio_Histo[i]->GetYaxis()->SetRangeUser(0,Ratio_p5_SM->GetBinContent(i)*2);
    ratio_Histo[i]->GetYaxis()->SetTitle("EFT/SM");
    if(EFT == "ctwi") ratio_Histo[i]->GetXaxis()->SetTitle("ctwi");
    if(EFT == "cbwi") ratio_Histo[i]->GetXaxis()->SetTitle("cbwi");
    ratio_Histo[i]->SetTitle("");
    Canvas->Print(file_name_eft.c_str());

    ratio_Histo[i]->Write(("bin_content_par1_"+number_plot).c_str());
    //ratio_formula->Write();
  }


  //file_output->Write();
  ratio_file->Close();




}

void Rwgt_vs(string variable, string EFT, int nbins, double xmin, double xmax, string selection, string W_value, TTree* t1, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string Name){

  //Fcts to plot base MC Simulations data vs Weighted MC simulation Data
  Name += "_" + variable; 

  TH1D *Histo_1 = new TH1D("Histo_1","Histo_1",nbins,xmin,xmax);
  TH1D* Histo_2 = GetHistoWeight(t1, variable, nbins, xmin, xmax, "1", "Histo_2");

  int nbfiles = 400;
  TFile* fInput[nbfiles];
  TTree* tInput[nbfiles];
  string inputName;

  string suffix[nbfiles];
  string file_name = "output_EVENTS_t_channel_madspin_Rwt_ctwi_cbwi_";

  //Alternative method
  /*for (int i=1; i<nbfiles; i++)
  {
    suffix[i] = to_string(i);
    inputName = "data/madgraph/weighted_events/root/output/" + file_name + suffix[i] + ".root";
    fInput[i] = new TFile(inputName.c_str(),"READ");
    tInput[i] = (TTree*) fInput[i]->Get("Tree");
    TH1D* Histo = GetHistoWeight(tInput[i],variable,nbins,xmin,xmax,selection,"");
    Histo_1->Add(Histo);
    Histo->Reset();
  }*/

  TList *list = new TList();
  for(int i=1 ; i<=nbfiles ; i++)
  {
    suffix[i] = to_string(i);
    inputName = "data/madgraph/weighted_events/root/output/" + file_name + suffix[i] + ".root"; 
    if(exists_test0(inputName))
    {
      cout<<inputName<<" found"<<endl;
      fInput[i] = new TFile(inputName.c_str(),"READ");
      tInput[i] = (TTree*) fInput[i]->Get("Tree");

      TH1D* Histo = GetHistoWeight(tInput[i],variable,nbins,xmin,xmax,selection,"");

      list->Add(Histo);
      Histo->Clear();  
    }
  }
  Histo_1->Reset();
  Histo_1->Merge(list);
  list->Clear();

 if(normalization)
 {
   double a = Histo_1->Integral();
   double b = Histo_2->Integral();
   Histo_1->Scale(1/a);
   Histo_2->Scale(1/b);

   Name = Name + "_normalized";
   legendY = legendY + "normalized";
 }

  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
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
   
   if(EFT=="ctwi")
   {
     legendEntry1 = "Rwgt C_{tW}^{I} = " + W_value + " TeV^{-2}";
     legendEntry2 = "C_{tW}^{I} = " + W_value + " TeV^{-2}";
   }
   
   if(EFT=="cbwi")
   {
     legendEntry1 = "Rwgt C_{bW}^{I} = " + W_value + " TeV^{-2}";
     legendEntry2 = "C_{bW}^{I} = " + W_value + " TeV^{-2}";
   }

    if(EFT == "SM")
   {
     legendEntry1 = "Rwgt SM";
     legendEntry2 = "SM";
   }



   TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
   legend->SetFillColor(kWhite);
   legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");
   legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(), "l");

   legend->Draw("SAME");

   cout << "Histo1 Entries: "<<Histo_1->Integral()<<endl;
   cout << "Histo2 Entries: "<<Histo_2->Integral()<<endl;

   Name = Name + "_" +selection + ".gif"; 

   Canvas->Print(Name.c_str());
}

void single_Rwgt(string variable, string EFT, int nbins, double xmin, double xmax, string selection, string W_value, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string Name){
  
  //Fcts to Plot Single weighted MC simulation Data
  Name += "_" + variable; 
  int nbfiles = 20;
  TH1D *Histo_1 = new TH1D("Histo_1","Histo_1",nbins,xmin,xmax);
  TFile* fInput[nbfiles];
  TTree* tInput[nbfiles];
  string inputName;

  string suffix[nbfiles];
  string file_name = "output_EVENTS_t_channel_madspin_Rwt_ctwi_cbwi_";

  //Alternative Way
  /*for (int i=1; i<=nbfiles; i++)
  {
    cout<<"Entered: "<<i<<endl;
    suffix[i] = to_string(i);
    inputName = "data/madgraph/weighted_events/root/output/" + file_name + suffix[i] + ".root";
    fInput[i] = new TFile(inputName.c_str(),"READ");
    tInput[i] = (TTree*) fInput[i]->Get("Tree");
    TH1D* Histo = GetHistoWeight(tInput[i],variable,nbins,xmin,xmax,selection,"");
    Histo_1->Add(Histo);
    Histo->Reset("ICES");
  }*/
 
  TList *list = new TList();
  for(int i=1 ; i<nbfiles ; i++)
  {
    suffix[i] = to_string(i);
    inputName = "data/madgraph/weighted_events/root/output/" + file_name + suffix[i] + ".root"; 
    if(exists_test0(inputName))
    {
      cout<<inputName<<" Exists!"<<endl;
      fInput[i] = new TFile(inputName.c_str(),"READ");
      tInput[i] = (TTree*) fInput[i]->Get("Tree");

      TH1D* Histo = GetHistoWeight(tInput[i],variable,nbins,xmin,xmax,selection,"");

      list->Add(Histo);
      Histo->Clear();
    }
  }
  Histo_1->Reset();
  Histo_1->Merge(list);
  list->Clear();
  cout << "Histo1 Entries: "<<Histo_1->Integral()<<endl;



 if(normalization)
 {
   double a = Histo_1->Integral();
   Histo_1->Scale(1/a);
   Name = Name + "_normalized";
   legendY = legendY + "normalized";
 } 

  double max = Histo_1->GetMaximum();
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
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

    if(EFT=="ctwi")
   {
     legendEntry1 = "Rwgt C_{tW}^{I} = " + W_value + " TeV^{-2}";
   }
   
   if(EFT=="cbwi")
   {
     legendEntry1 = "Rwgt C_{bW}^{I} = " + W_value + " TeV^{-2}";
   }

   if(EFT == "SM")
   {
     legendEntry1 = "Rwgt SM";
   }


   TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
   legend->SetFillColor(kWhite);
   legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");

   legend->Draw("SAME");


   Name = Name + "_" +selection + ".gif";
   Canvas->Print(Name.c_str());

}

void test_weights(string variable, string EFT, int nbins, double xmin, double xmax, string selection, string W_value, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string Name)
{
  //Fcts to Test Weights, Plot number of entries as a fcts of the nb of .root 
  //This fcts can only take up to 44 .rrot files [?], Memory issues?

  //TH1D *Histo_3 = new TH1D("histo","histo",200,0,200);

  Name = Name + "_" +selection;
  int nb_tests = 44;
  double mean = 0;
  double x[nb_tests] , y[nb_tests];
  for(int j=1 ; j<nb_tests ; j++)
  {
      int nbfiles = j;
      TH1D *Histo_1 = new TH1D("Histo_1","Histo_1",nbins,xmin,xmax);
      TFile* fInput[nbfiles];
      TTree* tInput[nbfiles];
      string inputName;

      string suffix[nbfiles];
      string file_name = "output_EVENTS_t_channel_madspin_Rwt_ctwi_cbwi_";



      TList *list = new TList();
      for(int i=1 ; i<nbfiles ; i++)
      {
        suffix[i] = to_string(i);
        inputName = "data/madgraph/weighted_events/root/output/" + file_name + suffix[i] + ".root"; 
        if(exists_test0(inputName))
        {
          //std::cout<<inputName << "Exists"<<std::endl;
          fInput[i] = new TFile(inputName.c_str(),"READ");
          tInput[i] = (TTree*) fInput[i]->Get("Tree");

          TH1D* Histo = GetHistoWeight(tInput[i],variable,nbins,xmin,xmax,selection,"");

          list->Add(Histo);
          Histo->Clear();
        }
  }
  Histo_1->Reset();
  Histo_1->Merge(list);
  list->Clear();

  //Histo_3->Fill(Histo_1->Integral());

  x[j]=Histo_1->Integral();
  mean+=x[j];
  y[j]=j;
  std::cout<<"x["<<j<<"]= "<<x[j]<<std::endl;

  Histo_1->Clear();

  }

  TCanvas *canvas1 = new TCanvas("c1","c1",600,600);

  TGraph *gr1 = new TGraph(nb_tests,x,y);
  gr1->GetXaxis()->SetLimits(x[0]-3,x[nb_tests-1]+3);
  gr1->GetYaxis()->SetLimits(0,100);
  gr1->SetMarkerStyle(kStar);
  gr1->SetTitle(selection.c_str());
  gr1->Draw("AP");
  //std::cout<<"Total value  = "<< mean <<std::endl;
  string meanV = "entries=" + to_string(mean);

  TLegend* legend = new TLegend(0.75, 0.75, 0.99, 0.99, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  //legend->AddEntry(gr1->GetName(), selection.c_str(), "l");
  legend->AddEntry(gr1->GetName(), meanV.c_str(), "l");

  legend->Draw("SAME");

  //Name+= "_" + meanV + ".gif";
  Name += ".gif";
  canvas1->SaveAs(Name.c_str());
}

void Compare_Weights(string variable, string EFT, int nbins, double xmin, double xmax, string selection, string W_value, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string Name){
  
  //Fcts to Plot all weights in the same Canvas

  //Name += "_" + variable; 
  string weights [9] = {"weight_ctwi_m2","weight_ctwi_m1","weight_ctwi_p1","weight_ctwi_p2","weight_cbwi_m2","weight_cbwi_m1","weight_cbwi_p1","weight_cbwi_p2","weight_SM"};

  int nbfiles = 201;
  TH1D *Histo_ctwi_m2 = new TH1D("ctwi_m2","ctwi_m2",nbins,xmin,xmax);
  TH1D *Histo_ctwi_m1 = new TH1D("ctwi_m1","ctwi_m1",nbins,xmin,xmax);
  TH1D *Histo_ctwi_p1 = new TH1D("ctwi_p1","ctwi_p1",nbins,xmin,xmax);
  TH1D *Histo_ctwi_p2 = new TH1D("ctwi_p2","ctwi_p2",nbins,xmin,xmax);
  TH1D *Histo_cbwi_m2 = new TH1D("cbwi_m2","cbwi_m2",nbins,xmin,xmax);
  TH1D *Histo_cbwi_m1 = new TH1D("cbwi_m1","cbwi_m1",nbins,xmin,xmax);
  TH1D *Histo_cbwi_p1 = new TH1D("cbwi_p1","cbwi_p1",nbins,xmin,xmax);
  TH1D *Histo_cbwi_p2 = new TH1D("cbwi_p2","cbwi_p2",nbins,xmin,xmax);
  TH1D *Histo_SM = new TH1D("SM","SM",nbins,xmin,xmax);


  TFile* fInput[nbfiles];
  TTree* tInput[nbfiles];
  string inputName;

  string suffix[nbfiles];
  string file_name = "output_EVENTS_t_channel_madspin_Rwt_ctwi_cbwi_";
 
  TList *list_ctwi_m2 = new TList();
  TList *list_ctwi_m1 = new TList();
  TList *list_ctwi_p1 = new TList();
  TList *list_ctwi_p2 = new TList();
  TList *list_cbwi_m2 = new TList();
  TList *list_cbwi_m1 = new TList();
  TList *list_cbwi_p1 = new TList();
  TList *list_cbwi_p2 = new TList();
  TList *list_SM = new TList();


  for(int i=1 ; i<nbfiles ; i++)
  {
    suffix[i] = to_string(i);
    inputName = "data/madgraph/weighted_events/root/output/" + file_name + suffix[i] + ".root"; 
    if(exists_test0(inputName))
    {
      //cout<<inputName<<" Exists!"<<endl;
      fInput[i] = new TFile(inputName.c_str(),"READ");
      tInput[i] = (TTree*) fInput[i]->Get("Tree");

      TH1D* Histo_ctwi_m2 = GetHistoWeight(tInput[i],"weight_ctwi_m2",nbins,xmin,xmax,"1","");
      TH1D* Histo_ctwi_m1 = GetHistoWeight(tInput[i],"weight_ctwi_m1",nbins,xmin,xmax,"1","");
      TH1D* Histo_ctwi_p1 = GetHistoWeight(tInput[i],"weight_ctwi_p1",nbins,xmin,xmax,"1","");
      TH1D* Histo_ctwi_p2 = GetHistoWeight(tInput[i],"weight_ctwi_p2",nbins,xmin,xmax,"1","");
      TH1D* Histo_cbwi_m2 = GetHistoWeight(tInput[i],"weight_cbwi_m2",nbins,xmin,xmax,"1","");
      TH1D* Histo_cbwi_m1 = GetHistoWeight(tInput[i],"weight_cbwi_m1",nbins,xmin,xmax,"1","");
      TH1D* Histo_cbwi_p1 = GetHistoWeight(tInput[i],"weight_cbwi_p1",nbins,xmin,xmax,"1","");
      TH1D* Histo_cbwi_p2 = GetHistoWeight(tInput[i],"weight_cbwi_p2",nbins,xmin,xmax,"1","");
      TH1D* Histo_SM = GetHistoWeight(tInput[i],"weight_SM",nbins,xmin,xmax,"1","");

      list_ctwi_m2->Add(Histo_ctwi_m2);
      list_ctwi_m1->Add(Histo_ctwi_m1);
      list_ctwi_p1->Add(Histo_ctwi_p1);
      list_ctwi_p2->Add(Histo_ctwi_p2);
      list_cbwi_m2->Add(Histo_cbwi_m2);
      list_cbwi_m1->Add(Histo_cbwi_m1);
      list_cbwi_p1->Add(Histo_cbwi_p1);
      list_cbwi_p2->Add(Histo_cbwi_p2);
      list_SM->Add(Histo_SM);
      
      Histo_ctwi_m2->Clear();
      Histo_ctwi_m1->Clear();
      Histo_ctwi_p1->Clear();
      Histo_ctwi_p2->Clear();
      Histo_cbwi_m2->Clear();
      Histo_cbwi_m1->Clear();
      Histo_cbwi_p1->Clear();
      Histo_cbwi_p2->Clear();
      Histo_SM->Clear();

    }
  }
  Histo_ctwi_m2->Reset();
  Histo_ctwi_m1->Reset();
  Histo_ctwi_p1->Reset();
  Histo_ctwi_p2->Reset();
  Histo_cbwi_m2->Reset();
  Histo_cbwi_m1->Reset();
  Histo_cbwi_p1->Reset();
  Histo_cbwi_p2->Reset();
  Histo_SM->Reset();

  Histo_ctwi_m2->Merge(list_ctwi_m2,"-NOL");
  Histo_ctwi_m1->Merge(list_ctwi_m1,"-NOL");
  Histo_ctwi_p1->Merge(list_ctwi_p1,"-NOL");
  Histo_ctwi_p2->Merge(list_ctwi_p2,"-NOL");
  Histo_cbwi_m2->Merge(list_cbwi_m2,"-NOL");
  Histo_cbwi_m1->Merge(list_ctwi_m1,"-NOL");
  Histo_cbwi_p1->Merge(list_cbwi_p1,"-NOL");
  Histo_cbwi_p2->Merge(list_cbwi_p2,"-NOL");
  Histo_SM->Merge(list_SM,"-NOL");

  double max = Histo_SM->GetMaximum();
  TCanvas* Canvas = new TCanvas("Canvas","Canvas",900,900);
  Histo_SM->SetTitle("");
  Histo_SM->SetAxisRange(0,max*2.2,"Y");
  //Histo_SM->SetAxisRange(-0.3,0.5,"X");
  //Histo_SM->SetAxisRange(-10,10,"X");
  Histo_SM->SetXTitle(legendX.c_str());
  Histo_SM->SetYTitle(legendY.c_str());
  Histo_SM->SetLineColor(kRed);
  Histo_SM->SetLineWidth(2);
  Histo_SM->Draw();

  Histo_ctwi_m2->SetLineColor(kBlue);
  Histo_ctwi_m2->SetLineWidth(2);
  Histo_ctwi_m2->Draw("SAME");

  Histo_ctwi_m1->SetLineColor(kGreen);
  Histo_ctwi_m1->SetLineWidth(2);
  Histo_ctwi_m1->Draw("SAME");

  Histo_ctwi_p1->SetLineColor(kViolet);
  Histo_ctwi_p1->SetLineWidth(2);
  Histo_ctwi_p1->Draw("SAME");

  Histo_ctwi_p2->SetLineColor(kYellow);
  Histo_ctwi_p2->SetLineWidth(2);
  Histo_ctwi_p2->Draw("SAME");

  Histo_cbwi_m2->SetLineColor(kBlue);
  Histo_cbwi_m2->SetLineStyle(3);
  Histo_cbwi_m2->SetLineWidth(2);
  Histo_cbwi_m2->Draw("SAME");

  Histo_cbwi_m1->SetLineColor(kGreen);
  Histo_cbwi_m1->SetLineStyle(3);
  Histo_cbwi_m1->SetLineWidth(2);
  Histo_cbwi_m1->Draw("SAME");

  Histo_cbwi_p1->SetLineColor(kViolet);
  Histo_cbwi_p1->SetLineStyle(3);
  Histo_cbwi_p1->SetLineWidth(2);
  Histo_cbwi_p1->Draw("SAME");

  Histo_cbwi_p2->SetLineColor(kYellow);
  Histo_cbwi_p2->SetLineStyle(3);
  Histo_cbwi_p2->SetLineWidth(2);
  Histo_cbwi_p2->Draw("SAME");


  double lx0, ly0, lx1, ly1;
   if (legendPlace=="legendUpLeft"){
 	 lx0 = 0.2;
 	 ly0 = 0.75;
 	 lx1 = 0.5;
 	 ly1 = 0.95;
   }
    if (legendPlace=="legendUpRight"){
 	 lx0 = 0.75;
 	 ly0 = 0.60;
 	 lx1 = 0.99;
 	 ly1 = 0.99;
   }


   TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
   legend->SetFillColor(kWhite);
   legend->AddEntry(Histo_SM->GetName(), "SM", "l");
   legend->AddEntry(Histo_ctwi_m2->GetName(), "C_{tw}^{I}=-2", "l");
   legend->AddEntry(Histo_ctwi_m1->GetName(), "C_{tw}^{I}=-1", "l");
   legend->AddEntry(Histo_ctwi_p1->GetName(), "C_{tw}^{I}=1", "l");
   legend->AddEntry(Histo_ctwi_p2->GetName(), "C_{tw}^{I}=2", "l");
   legend->AddEntry(Histo_cbwi_m2->GetName(), "C_{bw}^{I}=-2", "l"); 
   legend->AddEntry(Histo_cbwi_m1->GetName(), "C_{bw}^{I}=-1", "l"); 
   legend->AddEntry(Histo_cbwi_p1->GetName(), "C_{bw}^{I}=1", "l"); 
   legend->AddEntry(Histo_cbwi_p2->GetName(), "C_{bw}^{I}=2", "l"); 

   legend->Draw("SAME");


   Name = Name + "_" +selection + ".pdf";
   Canvas->Print(Name.c_str());

}


int main (){
int nbfiles = 13;
string suffix[nbfiles];
suffix[0] = "output_top_SM.root"; //SM

suffix[1] = "output_top_cbwi_n2.root"; // cbWi = -2
suffix[2] = "output_top_cbwi_n1.root";  //cbWi = -1
suffix[3] = "output_top_cbwi_p1.root"; //cbWi = 1
suffix[4] = "output_top_cbwi_p2.root"; //cbWi = 2

suffix[5] = "output_top_ctwi_n2.root"; //ctWi = -2
suffix[6] = "output_top_ctwi_n1.root"; //ctWi = -1
suffix[7] = "output_top_ctwi_p1.root"; //ctWi = 1
suffix[8] = "output_top_ctwi_p2.root"; //ctWi = 2

suffix[9] = "output_top_cbwi_n5.root"; //cbWi = -5
suffix[10] = "output_top_cbwi_p5.root"; //cbWi = 5
suffix[11] = "output_top_ctwi_n5.root"; //ctWi = -5
suffix[12] = "output_top_ctwi_p5.root"; //ctWi = 5

//suffix[13] = "output_madgraph_ctwi_m1.root";
//suffix[14] = "output_madgraph_ctwi_m2.root";
//suffix[15] = "output_madgraph_ctwi_p1.root";
//suffix[16] = "output_madgraph_ctwi_p2.root";
/*suffix[17] = "output_madspin_SM.root";
suffix[18] = "output_madspin_cbwi_m1.root";
suffix[19] = "output_madspin_cbwi_m2.root";
suffix[20] = "output_madspin_cbwi_p1.root";
suffix[21] = "output_madspin_cbwi_p2.root";
suffix[22] = "output_madspin_cptbi_m1.root";
suffix[23] = "output_madspin_cptbi_m2.root";
suffix[24] = "output_madspin_cptbi_p1.root";
suffix[25] = "output_madspin_cptbi_p.root";
suffix[26] = "output_madspin_ctw_m1.root";
suffix[27] = "output_madspin_ctw_m2.root";
suffix[28] = "output_madspin_ctw_p1.root";
suffix[29] = "output_madspin_ctw_p2.root";
suffix[30] = "output_madspin_ctwi_m1.root";
suffix[31] = "output_madspin_ctwi_m2.root";
suffix[32] = "output_madspin_ctwi_p1.root";
suffix[33] = "output_madspin_ctwi_p2.root";
suffix[34] = "/heppy/output_t_chan_MC.root";*/

TFile* fInput[nbfiles];
TTree* tInput[nbfiles];
string inputName;


for (int i=0; i<nbfiles; i++)
{
  inputName = "data/madgraph/output/1M/" + suffix[i];
  fInput[i] = new TFile(inputName.c_str(),"READ");
  tInput[i] = (TTree*) fInput[i]->Get("Tree");
}


  //////////cbWi Plots//////////
  //Compare_3Histos(tInput[0], tInput[1], tInput[4], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* [rad]", "", "legendUpRight", "dim6top", suffix[0], suffix[1], suffix[4], "cbWi", "results/cbWi_");
  //Compare_3Histos(tInput[0], tInput[1], tInput[4], "top_mass",20, 100, 400, "1", "Top Mass [GeV]", "", "legendUpRight", "dim6top", suffix[0], suffix[1], suffix[4], "cbWi", "results/cbWi_");
  //Compare_3Histos(tInput[0], tInput[1], tInput[4], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "", "legendUpRight", "dim6top", suffix[0], suffix[1], suffix[4], "cbWi", "results/cbWi_");

  //////////ctWi Plots//////////
  //Compare_3Histos(tInput[0], tInput[5], tInput[8], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* [rad]", "", "legendUpRight", "dim6top", suffix[0], suffix[5], suffix[8], "ctWi", "results/ctWi_");
  //Compare_3Histos(tInput[0], tInput[5], tInput[8], "top_mass",20, 100, 400, "1", "Top Mass [GeV]", "", "legendUpRight", "dim6top", suffix[0], suffix[5], suffix[8], "ctWi", "results/ctWi_");
  //Compare_3Histos(tInput[0], tInput[5], tInput[8], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "", "legendUpRight", "dim6top", suffix[0], suffix[5], suffix[8], "ctWi", "results/ctWi_");

  
  
if(EFT_vs_SM)
  {
    //costheta:
    //Ratio_EFT_SM(tInput[0],tInput[3],tInput[4],tInput[1],tInput[2],"cosTheta","cbwi",5,-1,1,"1","cos#theta","ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosTheta.pdf","lepton");
    
    //costhetaStar:
    //Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], "cosThetaStar","cbwi", 5, -1, 1 ,"nature_lepton==1", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosThetaStar.pdf", "elec");
    //Ratio_EFT_SM_7pts(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], tInput[10], tInput[9], "cosThetaStar","cbwi", 5, -1, 1 ,"nature_lepton==1", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosThetaStar_pt7.pdf", "elec");
    //Ratio_EFT_SM(tInput[0], tInput[7], tInput[8], tInput[5], tInput[6], "cosThetaStar","ctwi", 5, -1, 1 ,"nature_lepton==1", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_cosThetaStar.pdf", "elec");

    //PhiStar:
    //Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[2], tInput[1], "PhiStar","cbwi", 5, 0, 2*TMath::Pi(),"1", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_PhiStar.pdf","elec");
    //Ratio_EFT_SM(tInput[0], tInput[7], tInput[8], tInput[6], tInput[5], "PhiStar","ctwi", 5, 0, 2*TMath::Pi(),"1", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_PhiStar.pdf","elec");
    //Ratio_EFT_SM_7pts(tInput[0], tInput[7], tInput[8], tInput[6], tInput[5], tInput[12], tInput[11], "PhiStar","ctwi", 5, 0, 2*TMath::Pi(),"1", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_PhiStar.pdf","elec");
  } 

if(reweight)
  {
    /////////////////Change the variables to get the Plot you want/////////////////

    string EFT = "ctwi";                 //EFT variable {ctwi,cbwi,SM}
    string weight = "weight_ctwi_p2";   //Cut value
    int W_value = 2;                  //Wilson coefficient Value
    TTree *t1 = tInput[8];            //.root file to be analyzed

    Rwgt_vs("PhiStar", EFT, 20, 0, 6.5, weight, to_string(W_value), t1, "#phi^{*}", "", "legendUpRight", "dim6top", "", "", "results/weighted/Rwgt");
    Rwgt_vs("cosThetaStar", EFT, 20, -1, 1, weight, to_string(W_value), t1, "cos#theta*", "", "legendUpRight", "dim6top", "", "", "results/weighted/Rwgt");
    //Rwgt_vs("top_mass", EFT, 20, 168, 178, weight, to_string(W_value), t1, "Top Mass [GeV]", "", "legendUpRight", "dim6top", "", "", "results/weighted/Rwgt");
    //Rwgt_vs("lepton_pt", EFT, 20, 0, 100, weight, to_string(W_value), t1, "Lepton Pt [GeV]", "", "legendUpRight", "dim6top", "", "", "results/weighted/Rwgt");

    //single_Rwgt("PhiStar", EFT, 20, 0, 6.5, weight, to_string(W_value), "#phi^{*}", "", "legendUpRight", "dim6top", "", "results/weighted/single_Rwgt");
    
    //test_weights("PhiStar", EFT, 20, 0, 6.5, weight, to_string(W_value), "#phi^{*}", "", "legendUpRight", "dim6top", "", "results/weighted/test_weights/test_");
    
    //Compare_Weights("", EFT, 100, 0, 0.06, weight, to_string(W_value), "", "", "legendUpRight", "EFT Value [TeV^{-2}]", "", "results/weighted/compare_weightss");

  }

/*
  //Plot_SM
  //Compare_2Histos(tInput[0], tInput[17], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_CosTheta.pdf"  );
  Compare_2Histos(tInput[0], tInput[17], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_CosThetaStar.pdf");
  //Compare_2Histos(tInput[0], tInput[17], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_PhiStar.pdf");
  //Compare_2Histos(tInput[0], tInput[17], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_top_pt.pdf");
  Compare_2Histos(tInput[0], tInput[17], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_W_pt.pdf");
  Compare_2Histos(tInput[0], tInput[17], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[17], suffix[17], "results/dim6top_compareSM_lepton_pt.pdf");
	Compare_2Histos(tInput[0], tInput[17], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_lepton_E_Wframe.pdf");
  Compare_2Histos(tInput[0], tInput[17], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_top_mass.pdf");
  Compare_2Histos(tInput[0], tInput[17], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_W_mass.pdf");

  //MadGraph

  //Plots ctW
	Compare_3Histos(tInput[0], tInput[10], tInput[12], "cosTheta", 20, -1, 1, "1", "cos#theta", "a.u.", "legendUpLeft", "Operateur EFT", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_CosTheta.pdf");/*
	Compare_3Histos(tInput[0], tInput[10], tInput[12], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_CosThetaStar.pdf");
	Compare_3Histos(tInput[0], tInput[10], tInput[12], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_PhiStar.pdf");
  Compare_3Histos(tInput[0], tInput[10], tInput[12], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_top_pt.pdf");
  Compare_3Histos(tInput[0], tInput[10], tInput[12], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "Number of events", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_W_pt.pdf");
  /*Compare_3Histos(tInput[0], tInput[10], tInput[12], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_lepton_pt.pdf");
	Compare_3Histos(tInput[0], tInput[10], tInput[12], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[0], tInput[10], tInput[12], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "a.u.", "legendUpRight", "Operateur EFT", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_top_mass.pdf");
  /*Compare_3Histos(tInput[0], tInput[10], tInput[12], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_W_mass.pdf");
  Compare_3Histos(tInput[0], tInput[10], tInput[12], "W_transverse_mass",20, 50, 140, "1", "M_{T,W} (GeV)", "number of events", "legendUpRight", "Operateur de dimension 6", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_W_transverse_mass.pdf");

	//Plots ctWI*/
	//Compare_3Histos(tInput[0], tInput[14], tInput[16], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_CosTheta.pdf");
	//Compare_3Histos(tInput[0], tInput[14], tInput[16], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_CosThetaStar.pdf");
	//Compare_3Histos(tInput[0], tInput[14], tInput[16], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_PhiStar.pdf");
  /*(tInput[0], tInput[14], tInput[16], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_top_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[14], tInput[16], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_W_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[14], tInput[16], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_lepton_pt.pdf");
	//Compare_3Histos(tInput[0], tInput[14], tInput[16], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0],suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_lepton_E_Wframe.pdf");
  *///Compare_3Histos(tInput[0], tInput[14], tInput[16], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctwI_top_mass.pdf");
  /*Compare_3Histos(tInput[0], tInput[14], tInput[16], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_W_mass.pdf");

	//Plots cbWI
	//Compare_3Histos(tInput[0], tInput[2], tInput[4], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_CosTheta.pdf");
	Compare_3Histos(tInput[0], tInput[2], tInput[4], "cosThetaStar", 20, -1, 1, "1", "cos(#theta*)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_CosThetaStar.pdf");
	Compare_3Histos(tInput[0], tInput[2], tInput[4], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_PhiStar.pdf");
  Compare_3Histos(tInput[0], tInput[2], tInput[4], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_top_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[2], tInput[4], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_W_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[2], tInput[4], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_lepton_pt.pdf");
	//Compare_3Histos(tInput[0], tInput[2], tInput[4], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_lepton_E_Wframe.pdf");
  *///Compare_3Histos(tInput[0], tInput[2], tInput[4], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbwI_top_mass.pdf");
  /*Compare_3Histos(tInput[0], tInput[2], tInput[4], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbwI_W_mass.pdf");

  //Plots cptbI
	//Compare_3Histos(tInput[0], tInput[6], tInput[8], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_CosTheta.pdf");
	*///Compare_3Histos(tInput[0], tInput[6], tInput[8], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_CosThetaStar.pdf");
	//Compare_3Histos(tInput[0], tInput[6], tInput[8], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_PhiStar.pdf");
  //Compare_3Histos(tInput[0], tInput[6], tInput[8], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_top_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[6], tInput[8], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_W_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[6], tInput[8], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_lepton_pt.pdf");
	//Compare_3Histos(tInput[0], tInput[6], tInput[8], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_lepton_E_Wframe.pdf");
  //Compare_3Histos(tInput[0], tInput[6], tInput[8], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_top_mass.pdf");
  /*Compare_3Histos(tInput[0], tInput[6], tInput[8], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_W_mass.pdf");

  //MadSpin

  //Plots ctW
	Compare_3Histos(tInput[17], tInput[27], tInput[29], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_CosTheta.pdf");
	Compare_3Histos(tInput[17], tInput[27], tInput[29], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_CosThetaStar.pdf");
	Compare_3Histos(tInput[17], tInput[27], tInput[29], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_PhiStar.pdf");
  Compare_3Histos(tInput[17], tInput[27], tInput[29], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27],suffix[29], "results/madspin_dim6top_ctW_top_pt.pdf");
  Compare_3Histos(tInput[17], tInput[27], tInput[29], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_W_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[27], tInput[29], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_lepton_pt.pdf");
	//Compare_3Histos(tInput[17], tInput[27], tInput[29], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[17], tInput[27], tInput[29], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_top_mass.pdf");
  Compare_3Histos(tInput[17], tInput[27], tInput[29], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_W_mass.pdf");


  //Plots ctWI
	Compare_3Histos(tInput[17], tInput[31], tInput[33], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[31], suffix[33], "results/madspin_dim6top_ctWI_CosTheta.pdf");
	Compare_3Histos(tInput[17], tInput[31], tInput[33], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0],  suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_CosThetaStar.pdf");
	Compare_3Histos(tInput[17], tInput[31], tInput[33], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_PhiStar.pdf");
  Compare_3Histos(tInput[17], tInput[31], tInput[33], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0],  suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_top_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[31], tInput[33], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_W_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[31], tInput[33], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_lepton_pt.pdf");
	//Compare_3Histos(tInput[17], tInput[31], tInput[33], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0],  suffix[31], suffix[33], "results/madspin_dim6top_ctWI_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[17], tInput[31], tInput[33], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[31], suffix[33], "results/madspin_dim6top_ctWI_top_mass.pdf");
  Compare_3Histos(tInput[17], tInput[31], tInput[33], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[31], suffix[33], "results/madspin_dim6top_ctWI_W_mass.pdf");


	//Plots cbWI
	//Compare_3Histos(tInput[17], tInput[19], tInput[21], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_CosTheta.pdf");
	Compare_3Histos(tInput[17], tInput[19], tInput[21], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19],suffix[21], "results/madspin_dim6top_cbWI_CosThetaStar.pdf");
	Compare_3Histos(tInput[17], tInput[19], tInput[21], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_PhiStar.pdf");
  Compare_3Histos(tInput[17], tInput[19], tInput[21], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_top_pt.pdf");
  Compare_3Histos(tInput[17], tInput[19], tInput[21], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_W_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[19], tInput[21], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_lepton_pt.pdf");
	//Compare_3Histos(tInput[17], tInput[19], tInput[21], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[17], tInput[19], tInput[21], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbwI_top_mass.pdf");
  Compare_3Histos(tInput[17], tInput[19], tInput[21], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbwI_W_mass.pdf");

  //Plots cptbI
	//Compare_3Histos(tInput[17], tInput[23], tInput[25], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_CosTheta.pdf");
	//Compare_3Histos(tInput[17], tInput[23], tInput[25], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top",  suffix[0],  suffix[23], suffix[25], "results/madspin_dim6top_cptbI_CosThetaStar.pdf");
	Compare_3Histos(tInput[17], tInput[23], tInput[25], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_PhiStar.pdf");
  //Compare_3Histos(tInput[17], tInput[23], tInput[25], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_top_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[23], tInput[25], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[23],  suffix[25], "results/madspin_dim6top_cptbI_W_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[23], tInput[25], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_lepton_pt.pdf");
	//Compare_3Histos(tInput[17], tInput[23], tInput[25], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[17], tInput[23], tInput[25], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_top_mass.pdf");
  Compare_3Histos(tInput[17], tInput[23], tInput[25], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_W_mass.pdf");

  //Madgraph + MadSpin

  //cbwi
  //Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_CosTheta.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_CosThetaStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_PhiStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_top_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_W_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_lepton_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_lepton_E_Wframe.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_top_mass.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_W_mass.pdf");

  //cptbI
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_CosTheta.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_CosThetaStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_PhiStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_top_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_W_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_lepton_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_lepton_E_Wframe.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_top_mass.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_W_mass.pdf");

  //ctw
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_CosTheta.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_CosThetaStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_PhiStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_top_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_W_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_lepton_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_lepton_E_Wframe.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "Operateur de dimension 6", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_top_mass.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_W_mass.pdf");


  //ctwI
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_CosTheta.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_CosThetaStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_PhiStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_top_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_W_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_lepton_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_lepton_E_Wframe.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_top_mass.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_W_mass.pdf");
*/

  //-----------------Ratio EFT/SM-------------//

  //ctWI

  //Ratio_EFT_SM(tInput[0], tInput[15], tInput[16], tInput[13], tInput[14], "cosTheta","ctwi", 5, -1, 1,"1", "cos#theta", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_cosTheta.pdf");
  //Ratio_EFT_SM(tInput[0], tInput[15], tInput[16], tInput[13], tInput[14], "PhiStar","ctwi",20, 0, 6.2831 ,"nature_lepton == 1", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_PhiStar.pdf","elec");
  //Ratio_EFT_SM(tInput[0], tInput[15], tInput[16], tInput[13], tInput[14], "PhiStar","ctwi",20, 0, 6.2831 ,"nature_lepton == 2", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_PhiStar.pdf","muon");
  //Ratio_EFT_SM(tInput[0], tInput[15], tInput[16], tInput[13], tInput[14], "cosThetaStar","ctwi", 3, -1, 1 ,"1", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_cosThetaStar.pdf");

/*  Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], "cosTheta","cbwi", 5, -1, 1,"1", "cos#theta", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosTheta.pdf");
  Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], "PhiStar","cbwi", 5, 0, 2*TMath::Pi(),"1", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_PhiStar.pdf");*/
  //Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], "cosThetaStar","cbwi", 20, -1, 1 ,"nature_lepton==1", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosThetaStar.pdf", "elec");
  //Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], "cosThetaStar","cbwi", 20, -1, 1 ,"nature_lepton==2", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosThetaStar.pdf", "muon");
 
 
  /*Compare_3Histos(tInput[0], tInput[13], tInput[14], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[13], suffix[14], "results/madgraph_dim6top_ctWI_CosTheta_test.pdf");*/

  return 0;
}
