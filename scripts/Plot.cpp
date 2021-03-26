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
  
  //cout << "Area="<<histo->Integral()<<endl;
	//cout << "Nevents="<<histo->GetEntries()<<endl;
        histo->SetName(name.c_str());
        histo->SetTitle(name.c_str());

        return histo;
}

void Rwgt_vs(string selection, string variable, TTree* t1, TTree* t2,  int nbins, double xmin, double xmax, string normalized, string legendtitle, string legendX, string legendY, string legendPlace, string legendEntry1, string legendEntry2, string Name){

  //Fcts to plot base MC Simulations data vs Weighted MC simulation Data

  TH1D *Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D *Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, "1", "Histo_2");

 if(normalized == "normalized")
 {
   double a = Histo_1->Integral();
   double b = Histo_2->Integral();
   Histo_1->Scale(1/a);
   Histo_2->Scale(1/b);
 }
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");

  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();

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
   

   TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
   legend->SetFillColor(kWhite);

   legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");
   legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(), "l");

   legend->Draw("SAME");
   
   Canvas->Print(Name.c_str());
}

int main ()
{
int nbfiles = 14;
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

suffix[13] = "weighted_events.root"; //weighted Events

TFile* fInput[nbfiles];
TTree* tInput[nbfiles];
string inputName;


for (int i=0; i<nbfiles; i++)
{
  inputName = "data/madgraph/output/" + suffix[i];
  fInput[i] = new TFile(inputName.c_str(),"READ");
  tInput[i] = (TTree*) fInput[i]->Get("Tree");
}

//Phistar
Rwgt_vs("weight_SM", "PhiStar", tInput[13], tInput[0], 20, 0, 6.5, "normalized", "dim6top", "#phi^{*}", "normalized", "legendUpRight", "Rwgt SM", "SM", "results/Rwgt_Phistar_SM.gif");

Rwgt_vs("weight_cbwi_m2", "PhiStar", tInput[13], tInput[1], 20, 0, 6.5, "normalized", "dim6top", "#phi^{*}", "normalized", "legendUpRight", "Rwgt C_{bw}^{I} = -2", "C_{bw}^{I} = -2", "results/Rwgt_Phistar_cbwi_m2.gif");
Rwgt_vs("weight_cbwi_m1", "PhiStar", tInput[13], tInput[2], 20, 0, 6.5, "normalized", "dim6top", "#phi^{*}", "normalized", "legendUpRight", "Rwgt C_{bw}^{I} = -1", "C_{bw}^{I} = -1", "results/Rwgt_Phistar_cbwi_m1.gif");
Rwgt_vs("weight_cbwi_p1", "PhiStar", tInput[13], tInput[3], 20, 0, 6.5, "normalized", "dim6top", "#phi^{*}", "normalized", "legendUpRight", "Rwgt C_{bw}^{I} = 1", "C_{bw}^{I} = 1", "results/Rwgt_Phistar_cbwi_p1.gif");
Rwgt_vs("weight_cbwi_p2", "PhiStar", tInput[13], tInput[4], 20, 0, 6.5, "normalized", "dim6top", "#phi^{*}", "normalized", "legendUpRight", "Rwgt C_{bw}^{I} = 2", "C_{bw}^{I} = 2", "results/Rwgt_Phistar_cbwi_p2.gif");

Rwgt_vs("weight_ctwi_m2", "PhiStar", tInput[13], tInput[5], 20, 0, 6.5, "normalized", "dim6top", "#phi^{*}", "normalized", "legendUpRight", "Rwgt C_{tw}^{I} = -2", "C_{tw}^{I} = -2", "results/Rwgt_Phistar_ctwi_m2.gif");
Rwgt_vs("weight_ctwi_m1", "PhiStar", tInput[13], tInput[6], 20, 0, 6.5, "normalized", "dim6top", "#phi^{*}", "normalized", "legendUpRight", "Rwgt C_{tw}^{I} = -1", "C_{tw}^{I} = -1", "results/Rwgt_Phistar_ctwi_m1.gif");
Rwgt_vs("weight_ctwi_p1", "PhiStar", tInput[13], tInput[7], 20, 0, 6.5, "normalized", "dim6top", "#phi^{*}", "normalized", "legendUpRight", "Rwgt C_{tw}^{I} = 1", "C_{tw}^{I} = 1", "results/Rwgt_Phistar_ctwi_p1.gif");
Rwgt_vs("weight_ctwi_p2", "PhiStar", tInput[13], tInput[8], 20, 0, 6.5, "normalized", "dim6top", "#phi^{*}", "normalized", "legendUpRight", "Rwgt C_{tw}^{I} = 2", "C_{tw}^{I} = 2", "results/Rwgt_Phistar_ctwi_p2.gif");


//cosThetaStar
Rwgt_vs("weight_cbwi_m2", "cosThetaStar", tInput[13], tInput[1], 20, -1, 1, "normalized", "dim6top", "cos(#theta^{*})", "normalized", "legendUpRight", "Rwgt C_{bw}^{I} = -2", "C_{bw}^{I} = -2", "results/Rwgt_cosThetaStar_cbwi_m2.gif");

//Top Mass
Rwgt_vs("weight_cbwi_m2", "top_mass", tInput[13], tInput[1], 20, 100, 400, "normalized", "dim6top", "Top mass [GeV]", "normalized", "legendUpRight", "Rwgt C_{bw}^{I} = -2", "C_{bw}^{I} = -2", "results/Rwgt_TopMass_cbwi_m2.gif");



return 0;

}
