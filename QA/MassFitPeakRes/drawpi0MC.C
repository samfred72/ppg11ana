#include <iostream>
#include <RooWorkspace.h>
#include <cmath>
#include "../../histmaking/ana.cxx"
#include "../../histmaking/TreeSetting.h"
#include "/sphenix/user/jpark4/Utility/Style_jaebeom.h"
#include "/sphenix/user/jpark4/Utility/commonUtility.h"
//#include "getEff.cc"

  
double GetGraphMaximum(TGraph* graph) {
  return graph ? TMath::MaxElement(graph->GetN(), graph->GetY()) : 0;
}

void drawpi0MC()
{
  std::string particle="pi0";
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);
  ana anac;
  const int nbins = anac.nPtBins;
  const double* bins = anac.ptBins;
  double sPHENIX_posx = anac.sPHENIX_posx;
  double sPHENIX_posy = anac.sPHENIX_posy;
  double posy_diff = anac.posy_diff;

  TGraphErrors* gmc_mean = new TGraphErrors();
  TGraphErrors* gmc_res = new TGraphErrors();
  TGraphErrors* gmc_ratio = new TGraphErrors();
  
  int pointcount=0;
  for(int i = 1; i < 7; ++i) {
    TFile *f = new TFile(Form("/sphenix/user/samfred/run25/ppg11/histmaking/sufficient/workspace_fits_MB_%s_pt%i.root",particle.c_str(),i),"read");

    RooWorkspace *w = (RooWorkspace*) f->Get("workspace");
    RooFitResult *fitRes = (RooFitResult*)w->obj("fitresult_model_binnedData");

    //MC
    RooRealVar *mean = w->var("mean"); 
    RooRealVar *sigma = w->var("sigma1"); 
    RooFormulaVar * ratio = new RooFormulaVar("ratio", "sigma1 / mean", RooArgList(*sigma, *mean));

    gmc_mean->SetPoint(pointcount, (i+0.5), mean->getVal());
    gmc_mean->SetPointError(pointcount, 0.5, mean->getError());
    
    gmc_res->SetPoint(pointcount, (i+0.5), sigma->getVal());
    gmc_res->SetPointError(pointcount, 0.5, sigma->getError());
    
    gmc_ratio->SetPoint(pointcount, (i+0.5), ratio->getVal());
    gmc_ratio->SetPointError(pointcount, 0.5, ratio->getPropagatedError(*fitRes));

    f->Close();
    
    pointcount++;
  }
  
  gmc_mean->GetXaxis()->SetLimits(0,10);
  gmc_mean->GetXaxis()->SetRangeUser(0,10);
  gmc_mean->GetXaxis()->SetLimits(0,10);
  gmc_mean->GetXaxis()->SetRangeUser(0,10);
  gmc_res->GetXaxis()->SetLimits(0,10);
  gmc_res->GetXaxis()->SetRangeUser(0,10);
  gmc_res->GetXaxis()->SetLimits(0,10);
  gmc_res->GetXaxis()->SetRangeUser(0,10);
  gmc_ratio->GetXaxis()->SetLimits(0,10);
  gmc_ratio->GetXaxis()->SetRangeUser(0,10);
  gmc_ratio->GetXaxis()->SetLimits(0,10);
  gmc_ratio->GetXaxis()->SetRangeUser(0,10);

  SetGraphStyle(gmc_mean,0,0);
  SetGraphStyle(gmc_res,0,0);
  SetGraphStyle(gmc_ratio,0,0);

  gmc_mean->SetMarkerStyle(kOpenCircle);
  gmc_res->SetMarkerStyle(kOpenCircle);
  gmc_ratio->SetMarkerStyle(kOpenCircle);
  
  gmc_mean->GetYaxis()->SetLimits(0.05,0.3);
  gmc_mean->GetYaxis()->SetRangeUser(0.11,0.2);
  gmc_res->GetYaxis()->SetLimits(0.0,0.05);
  gmc_res->GetYaxis()->SetRangeUser(0.0,0.05);
  gmc_ratio->GetYaxis()->SetLimits(0.0,0.5);
  gmc_ratio->GetYaxis()->SetRangeUser(0.08,0.2);


  // mean
  TCanvas* c1 = new TCanvas("c1","",700,700);
  c1->cd();
  c1->SetTicks(1,1);
  c1->SetRightMargin(0.10);
  c1->SetLeftMargin(0.13);
  c1->SetTopMargin(0.06);
  gmc_mean->Draw("AP");
  gmc_mean->GetXaxis()->SetTitle("p_{T}^{#pi^{0}} [GeV]");
  gmc_mean->GetYaxis()->SetTitle("m_{0} [GeV]");
  gmc_mean->GetYaxis()->CenterTitle();
  gmc_mean->GetXaxis()->CenterTitle();
  
  drawText("#bf{#it{sPHENIX}} Internal",sPHENIX_posx,sPHENIX_posy,1,22);
  drawText("pp #sqrt{s} = 200 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff,1,19);
  drawText("|v_{z}^{MBD}| < 30 cm",sPHENIX_posx,sPHENIX_posy-posy_diff*2,1,18);
  drawText("p_{T}^{#gamma} > 0.5 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff*3,1,18);
  drawText("#alpha < 0.6",sPHENIX_posx,sPHENIX_posy-posy_diff*4,1,18);
  drawText("#gamma prob. > 0.05",sPHENIX_posx,sPHENIX_posy-posy_diff*5,1,18);
  drawText("#pi^{0}",sPHENIX_posx,sPHENIX_posy-posy_diff*6,1,18);

  TLegend *l1 = new TLegend(0.18,0.71,0.37,0.88);
  SetLegendStyle(l1);
  l1->SetTextSize(0.027);
  l1->SetHeader("Pythia Run 28");
  l1->AddEntry(gmc_mean,"Pythia MB, p_{T}^{leading #gamma} > 0.5 GeV","pe");
  l1->Draw("same");
  
  c1->SaveAs("plot_pi0_peak_mc.pdf");

  // width
  TCanvas* c2 = new TCanvas("c2","",700,700);
  c2->cd();
  c2->SetTicks(1,1);
  c2->SetRightMargin(0.10);
  c2->SetLeftMargin(0.13);
  c2->SetTopMargin(0.06);
  gmc_res->Draw("AP");
  gmc_res->GetXaxis()->SetTitle("p_{T}^{#pi^{0}} [GeV]");
  gmc_res->GetYaxis()->SetTitle("m_{0} width [GeV]");
  gmc_res->GetYaxis()->CenterTitle();
  gmc_res->GetXaxis()->CenterTitle();
  
  drawText("#bf{#it{sPHENIX}} Internal",sPHENIX_posx,sPHENIX_posy,1,22);
  drawText("pp #sqrt{s} = 200 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff,1,19);
  drawText("|v_{z}^{MBD}| < 30 cm",sPHENIX_posx,sPHENIX_posy-posy_diff*2,1,18);
  drawText("p_{T}^{#gamma} > 0.5 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff*3,1,18);
  drawText("#alpha < 0.6",sPHENIX_posx,sPHENIX_posy-posy_diff*4,1,18);
  drawText("#gamma prob. > 0.05",sPHENIX_posx,sPHENIX_posy-posy_diff*5,1,18);
  drawText("#pi^{0}",sPHENIX_posx,sPHENIX_posy-posy_diff*6,1,18);
  
  TLegend *l2 = new TLegend(0.18,0.71,0.37,0.88);
  SetLegendStyle(l2);
  l2->SetTextSize(0.027);
  l2->SetHeader("Pythia Run 28");
  l2->AddEntry(gmc_res,"Pythia MB, p_{T}^{leading #gamma} > 0.5 GeV","pe");
  l2->Draw("same");
  

  // ratio
  TCanvas* c3 = new TCanvas("c3","",700,700);
  c3->cd();
  c3->SetTicks(1,1);
  c3->SetRightMargin(0.10);
  c3->SetLeftMargin(0.13);
  c3->SetTopMargin(0.06);
  gmc_ratio->Draw("AP");
  gmc_ratio->GetXaxis()->SetTitle("p_{T}^{#pi^{0}} [GeV]");
  gmc_ratio->GetYaxis()->SetTitle("width/m_{0}");
  gmc_ratio->GetYaxis()->CenterTitle();
  gmc_ratio->GetXaxis()->CenterTitle();
  
  drawText("#bf{#it{sPHENIX}} Internal",sPHENIX_posx,sPHENIX_posy,1,22);
  drawText("pp #sqrt{s} = 200 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff,1,19);
  drawText("|v_{z}^{MBD}| < 30 cm",sPHENIX_posx,sPHENIX_posy-posy_diff*2,1,18);
  drawText("p_{T}^{#gamma} > 0.5 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff*3,1,18);
  drawText("#alpha < 0.6",sPHENIX_posx,sPHENIX_posy-posy_diff*4,1,18);
  drawText("#gamma prob. > 0.05",sPHENIX_posx,sPHENIX_posy-posy_diff*5,1,18);
  drawText("#pi^{0}",sPHENIX_posx,sPHENIX_posy-posy_diff*6,1,18);
  
  TLegend *l3 = new TLegend(0.18,0.71,0.37,0.88);
  SetLegendStyle(l3);
  l3->SetTextSize(0.027);
  l3->SetHeader("Pythia Run 28");
  l3->AddEntry(gmc_res,"Pythia MB, p_{T}^{leading #gamma} > 0.5 GeV","pe");
  l3->Draw("same");

}
