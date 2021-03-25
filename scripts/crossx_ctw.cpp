#include <iostream>
#include <TH2F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
using namespace std;

int main ()
{
  TCanvas* support = new TCanvas("","graph",450,450);
  support->Draw();

  TGraph *axe = new TGraph (5);
  TGraph *axe1 = new TGraph (5);
  TGraph *axe2 = new TGraph (5);
  TGraph *axe3 = new TGraph (5);
  TGraph *axe4 = new TGraph (5);
  TGraph *axe5 = new TGraph (5);
  TGraph *axe6 = new TGraph (5);
  TGraph *axe7 = new TGraph (5);


  auto legend = new TLegend(0.70,0.20,0.95,0.45);
  TH2F* axis = new TH2F("","",10,-3,3,10,30,50);
  //axis->SetTitle("Cross section as function of tw coupling with complex mass scheme at #delta#sigma #propto 10^{-2} ");
  axis->GetYaxis()->SetTitle("#sigma (pb)");
  axis->GetXaxis()->SetTitle("c_{tw}");
  axis->SetStats(kFALSE);
  axis->Draw();


  //Halfchain_ctw

  axe->SetLineColor(kBlue);
  axe->SetPoint(0,-2,41.1);
  axe->SetPoint(1,-1,40.4);
  axe->SetPoint(2,0,40.79);
  axe->SetPoint(3,1,42.29);
  axe->SetPoint(4,2,44.92);
  axe->Fit("pol2","","",-3,3);
  axe->GetFunction("pol2")->SetLineColor(kBlue);
  axe->Draw("*SAME");

  //Halfchain_ctwI

  axe5->SetPoint(0,-2,43.25);
  axe5->SetPoint(1,-1,41.45);
  axe5->SetPoint(2,0,40.79);
  axe5->SetPoint(3,1,41.25);
  axe5->SetPoint(4,2,42.84);
  axe5->Fit("pol2","","",-3,3);
  axe5->GetFunction("pol2")->SetLineColor(kBlack);
  //axe5->Draw("*SAME");

  //Halfchain_cbwI

  axe6->SetPoint(0,-2,43.06);
  axe6->SetPoint(1,-1,41.36);
  axe6->SetPoint(2,0,40.79);
  axe6->SetPoint(3,1,41.36);
  axe6->SetPoint(4,2,43.06);
  axe6->Fit("pol2","","",-3,3);
  axe6->GetFunction("pol2")->SetLineColor(kViolet);
  //axe6->Draw("*SAME");

  //Halfchain_cptbI

  axe7->SetPoint(0,-2,40.92);
  axe7->SetPoint(1,-1,40.83);
  axe7->SetPoint(2,0,40.79);
  axe7->SetPoint(3,1,40.83);
  axe7->SetPoint(4,2,40.92);
  axe7->Fit("pol2","","",-3,3);
  axe7->GetFunction("pol2")->SetLineColor(kGreen);
  //axe7->Draw("*SAME");

  //Fullchain_cbwI

  axe1->SetPoint(0,-2,42.59);
  axe1->SetPoint(1,-1,40.9);
  axe1->SetPoint(2,0,40.32);
  axe1->SetPoint(3,1,40.9);
  axe1->SetPoint(4,2,42.59);
  axe1->Fit("pol2","","",-3,3);
  axe1->GetFunction("pol2")->SetLineColor(kBlack);
  //axe1->Draw("*SAME");

  //Fullchain_ctw

  axe2->SetPoint(0,-2,40.71);
  axe2->SetPoint(1,-1,39.98);
  axe2->SetPoint(2,0,40.32);
  axe2->SetPoint(3,1,41.82);
  axe2->SetPoint(4,2,44.39);
  axe2->Fit("pol2","","",-3,3);
  axe2->GetFunction("pol2")->SetLineColor(kRed);
  axe2->Draw("*SAME");

  //Fullchain_ctwI

  axe3->SetPoint(0,-2,42.47);
  axe3->SetPoint(1,-1,40.99);
  axe3->SetPoint(2,0,40.32);
  axe3->SetPoint(3,1,40.8);
  axe3->SetPoint(4,2,42.37);
  axe3->Fit("pol2","","",-3,3);
  axe3->GetFunction("pol2")->SetLineColor(kOrange);
  //axe3->Draw("*SAME");

  //Fullchain_cptbI

  axe4->SetPoint(0,-2,40.47);
  axe4->SetPoint(1,-1,40.37);
  axe4->SetPoint(2,0,40.32);
  axe4->SetPoint(3,1,40.37);
  axe4->SetPoint(4,2,40.47);
  axe4->Fit("pol2","","",-3,3);
  axe4->GetFunction("pol2")->SetLineColor(kViolet);
  //axe4->Draw("*SAME");


  legend->AddEntry(axe->GetFunction("pol2"), "madspin","l");
  //legend->AddEntry(axe1->GetFunction("pol2"), "madgraph","l");
  legend->AddEntry(axe2->GetFunction("pol2"), "madgraph","l");
  //legend->AddEntry(axe3->GetFunction("pol2"), "madgraph","l");
  //legend->AddEntry(axe4->GetFunction("pol2"), "madgraph","l");
  //legend->AddEntry(axe5->GetFunction("pol2"), "madspin","l");
  //legend->AddEntry(axe6->GetFunction("pol2"), "madspin","l");
  //legend->AddEntry(axe7->GetFunction("pol2"), "madspin","l");

  legend->Draw("SAME");
  support->SaveAs("results/ctw_comparision.pdf");
  return 0;
}
