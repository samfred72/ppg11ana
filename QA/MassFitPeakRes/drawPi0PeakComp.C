#include <iostream>
#include <RooWorkspace.h>
#include <cmath>
#include "../Headers/ana.cxx"
#include "../Headers/TreeSetting.h"
#include "/sphenix/user/jpark4/Utility/Style_jaebeom.h"
#include "/sphenix/user/jpark4/Utility/commonUtility.h"
//#include "getEff.cc"

  
double GetGraphMaximum(TGraph* graph) {
  return graph ? TMath::MaxElement(graph->GetN(), graph->GetY()) : 0;
}

void drawPi0PeakComp()
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

  TGraphErrors* gdata_mean1 = new TGraphErrors();
  TGraphErrors* gdata_mean2 = new TGraphErrors();
  TGraphErrors* gmc_mean1 = new TGraphErrors();
  TGraphErrors* gmc_mean2 = new TGraphErrors();
  TGraphErrors* gmc_mean3 = new TGraphErrors();

  TGraphErrors* gdata_res1 = new TGraphErrors();
  TGraphErrors* gdata_res2 = new TGraphErrors();
  TGraphErrors* gmc_res1 = new TGraphErrors();
  TGraphErrors* gmc_res2 = new TGraphErrors();
  
  int pointcount=0;
  for(int i=0; i < 9; ++i){
    float ptlow = bins[i];
    float pthigh = bins[i+1];
    TFile *f1 = new TFile(Form("../Fitting/fitresultsGC/fitresGC_pi0_DoubleGaus_cheb3_split_trigopt2_ptLow%.1f_ptHigh%.1f.root",ptlow,pthigh),"read");
    TFile *f2 = new TFile(Form("../Fitting/fitresultsMC/fitresMC_pi0_DoubleGaus_split_trigopt2_ptLow%.1f_ptHigh%.1f.root",ptlow,pthigh),"read");
    TFile *f3 = new TFile(Form("../Fitting/fitresultsMC/weighted_fitresMC_pi0_DoubleGaus_split_trigopt2_ptLow%.1f_ptHigh%.1f.root",ptlow,pthigh),"read");

    RooWorkspace *w1 = (RooWorkspace*) f1->Get("workspace");
    RooWorkspace *w2 = (RooWorkspace*) f2->Get("workspace");
    RooWorkspace *w3 = (RooWorkspace*) f3->Get("workspace");
    RooFitResult *fitResData = (RooFitResult*)w1->obj("fitresult_modelgc_binnedData");
    RooFitResult *fitResMC = (RooFitResult*)w2->obj("fitresult_model_binnedData");
    
    //DATA
    RooRealVar *meandata = w1->var("mean"); 
    RooRealVar *sigmadata = w1->var("sigma1"); 
    RooRealVar *xdata = w1->var("x"); 
    RooRealVar *fracdata = w1->var("frac"); 
    meandata->SetName("meandata");
    sigmadata->SetName("sigmadata");
    xdata->SetName("xdata");
    fracdata->SetName("fracdata");

    RooFormulaVar sigmaresdata("sigmaresdata", "(sigmadata*fracdata + sigmadata*xdata*(1-fracdata)) / meandata", RooArgList(*sigmadata, *fracdata, *xdata, *meandata));
    
    double sigmaresval = sigmaresdata.getVal();
    double sigmareserr = sigmaresdata.getPropagatedError(*fitResData);

    gdata_mean2->SetPoint(pointcount, (ptlow+pthigh)/2, meandata->getVal());
    gdata_mean2->SetPointError(pointcount, (pthigh-ptlow)/2, meandata->getError());
    
    gdata_res2->SetPoint(pointcount, (ptlow+pthigh)/2, sigmaresval); 
    gdata_res2->SetPointError(pointcount, (pthigh-ptlow)/2, sigmareserr);

    //MC
    RooRealVar *meanmc = w2->var("mean"); 
    RooRealVar *sigmamc = w2->var("sigma1"); 
    RooRealVar *xmc = w2->var("x"); 
    RooRealVar *fracmc = w2->var("frac"); 
    RooFormulaVar sigmaresmc("sigmaresmc", "(sigma1*frac + sigma1*x*(1-frac)) / mean", RooArgList(*sigmamc, *fracmc, *xmc, *meanmc));
    
    sigmaresval = sigmaresmc.getVal();
    sigmareserr = sigmaresmc.getPropagatedError(*fitResMC);

    if(ptlow==5){
      std::cout << sigmamc->getVal() << " +/- " << sigmamc->getError() << std::endl;
      std::cout << xmc->getVal() << " +/- " << xmc->getError() << std::endl;
      std::cout << fracmc->getVal() << " +/- " << fracmc->getError() << std::endl;
      std::cout << sigmareserr << std::endl;
    } 

    gmc_mean2->SetPoint(pointcount, (ptlow+pthigh)/2, meanmc->getVal());
    gmc_mean2->SetPointError(pointcount, (pthigh-ptlow)/2, meanmc->getError());
    
    RooRealVar *meanmcw = w3->var("mean"); 
    gmc_mean3->SetPoint(pointcount, (ptlow+pthigh)/2, meanmcw->getVal());
    gmc_mean3->SetPointError(pointcount, (pthigh-ptlow)/2, 0);
    gmc_res2->SetPoint(pointcount, (ptlow+pthigh)/2, sigmaresval); 
    gmc_res2->SetPointError(pointcount, (pthigh-ptlow)/2, sigmareserr);
    
    f1->Close();
    f2->Close();
    
    pointcount++;
  }
 

  pointcount=0;
  for(int i=8; i < 12; ++i){
    float ptlow = bins[i];
    float pthigh = bins[i+1];
    TFile *f1 = new TFile(Form("../Fitting/fitresultsGC/fitresGC_pi0_DoubleGaus_cheb3_split_trigopt1_ptLow%.1f_ptHigh%.1f.root",ptlow,pthigh),"read");
    TFile *f2 = new TFile(Form("../Fitting/fitresultsMC/fitresMC_pi0_DoubleGaus_split_trigopt1_ptLow%.1f_ptHigh%.1f.root",ptlow,pthigh),"read");

    if(!f1 || !f2) {std::cout << "error..." << std::endl;}

    RooWorkspace *w1 = (RooWorkspace*) f1->Get("workspace");
    RooWorkspace *w2 = (RooWorkspace*) f2->Get("workspace");
    RooFitResult *fitResData = (RooFitResult*)w1->obj("fitresult_modelgc_binnedData");
    RooFitResult *fitResMC = (RooFitResult*)w2->obj("fitresult_model_binnedData");
    
    //DATA
    RooRealVar *meandata = w1->var("mean"); 
    RooRealVar *sigmadata = w1->var("sigma1"); 
    RooRealVar *xdata = w1->var("x"); 
    RooRealVar *fracdata = w1->var("frac"); 
    meandata->SetName("meandata");
    sigmadata->SetName("sigmadata");
    xdata->SetName("xdata");
    fracdata->SetName("fracdata");

    RooFormulaVar sigmaresdata("sigmaresdata", "(sigmadata*fracdata + sigmadata*xdata*(1-fracdata)) / meandata", RooArgList(*sigmadata, *fracdata, *xdata, *meandata));
    

    double sigmaresval = sigmaresdata.getVal();
    double sigmareserr = sigmaresdata.getPropagatedError(*fitResData);

    gdata_mean1->SetPoint(pointcount, (ptlow+pthigh)/2, meandata->getVal());
    gdata_mean1->SetPointError(pointcount, (pthigh-ptlow)/2, meandata->getError());
    
    gdata_res1->SetPoint(pointcount, (ptlow+pthigh)/2, sigmaresval); 
    gdata_res1->SetPointError(pointcount, (pthigh-ptlow)/2, sigmareserr);

    //MC
    RooRealVar *meanmc = w2->var("mean"); 
    RooRealVar *sigmamc = w2->var("sigma1"); 
    RooRealVar *xmc = w2->var("x"); 
    RooRealVar *fracmc = w2->var("frac"); 

    RooFormulaVar sigmaresmc("sigmaresmc", "(sigma1*frac + sigma1*x*(1-frac)) / mean", RooArgList(*sigmamc, *fracmc, *xmc, *meanmc));

    sigmaresval = sigmaresmc.getVal();
    sigmareserr = sigmaresmc.getPropagatedError(*fitResMC);

    gmc_mean1->SetPoint(pointcount, (ptlow+pthigh)/2, meanmc->getVal());
    gmc_mean1->SetPointError(pointcount, (pthigh-ptlow)/2, meanmc->getError());
    
    gmc_res1->SetPoint(pointcount, (ptlow+pthigh)/2, sigmaresval); 
    gmc_res1->SetPointError(pointcount, (pthigh-ptlow)/2, sigmareserr);
    f1->Close();
    f2->Close();
    pointcount++;
  }
  
  gdata_mean1->GetXaxis()->SetLimits(0,10);
  gdata_mean1->GetXaxis()->SetRangeUser(0,10);
  gdata_mean1->GetXaxis()->SetLimits(0,10);
  gdata_mean1->GetXaxis()->SetRangeUser(0,10);

  gdata_res1->GetXaxis()->SetLimits(0,10);
  gdata_res1->GetXaxis()->SetRangeUser(0,10);
  gdata_res1->GetXaxis()->SetLimits(0,10);
  gdata_res1->GetXaxis()->SetRangeUser(0,10);

  SetGraphStyle(gdata_mean1,0,0);
  SetGraphStyle(gdata_mean2,0,0);

  SetGraphStyle(gmc_mean1,1,1);
  SetGraphStyle(gmc_mean2,1,1);
  SetGraphStyle(gmc_mean3,2,1);

  SetGraphStyle(gdata_res1,0,0);
  SetGraphStyle(gdata_res2,0,0);

  SetGraphStyle(gmc_res1,1,1);
  SetGraphStyle(gmc_res2,1,1);

  gdata_mean2->SetMarkerStyle(kOpenCircle);
  gmc_mean2->SetMarkerStyle(kOpenCircle);
  gmc_mean3->SetMarkerStyle(kOpenCircle);
  gdata_res2->SetMarkerStyle(kOpenCircle);
  gmc_res2->SetMarkerStyle(kOpenCircle);
  
  gdata_mean1->GetYaxis()->SetLimits(0.05,0.3);
  gdata_mean1->GetYaxis()->SetRangeUser(0.11,0.2);
  
  gdata_res1->GetYaxis()->SetLimits(0.0,0.4);
  gdata_res1->GetYaxis()->SetRangeUser(0.0,0.4);


  TCanvas* c1 = new TCanvas("c1","",700,700);
  c1->cd();
  c1->SetTicks(1,1);
  c1->SetRightMargin(0.10);
  c1->SetLeftMargin(0.13);
  c1->SetTopMargin(0.06);
  gdata_mean1->Draw("AP");
  gdata_mean2->Draw("P same");
  gmc_mean1->Draw("P same");
  gmc_mean2->Draw("P same");
  gmc_mean3->Draw("P same");
  gdata_mean1->GetXaxis()->SetTitle("p_{T}^{#pi^{0}} [GeV]");
  gdata_mean1->GetYaxis()->SetTitle("m_{0} [GeV]");
  gdata_mean1->GetYaxis()->CenterTitle();
  gdata_mean1->GetXaxis()->CenterTitle();
  
  drawText("#bf{#it{sPHENIX}} Internal",sPHENIX_posx,sPHENIX_posy,1,22);
  drawText("pp #sqrt{s} = 200 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff,1,19);
  drawText("|v_{z}^{MBD}| < 30 cm",sPHENIX_posx,sPHENIX_posy-posy_diff*2,1,18);
  drawText("p_{T}^{#gamma} > 1 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff*3,1,18);
  drawText("#alpha < 0.6",sPHENIX_posx,sPHENIX_posy-posy_diff*4,1,18);
  drawText("#gamma prob. > 0.05",sPHENIX_posx,sPHENIX_posy-posy_diff*5,1,18);
  drawText("#pi^{0}",sPHENIX_posx,sPHENIX_posy-posy_diff*6,1,18);

  TLegend *l1 = new TLegend(0.18,0.71,0.37,0.88);
  SetLegendStyle(l1);
  l1->SetTextSize(0.027);
  l1->SetHeader("Data");
  l1->AddEntry(gdata_mean1,"Photon 4 GeV, p_{T}^{leading #gamma} > 6 GeV","pe");
  l1->AddEntry(gdata_mean2,"MBD N&S>=1, p_{T}^{leading #gamma} > 1 GeV","pe");
  l1->Draw("same");
  
  TLegend *l2 = new TLegend(0.18,0.52,0.37,0.69);
  SetLegendStyle(l2);
  l2->SetTextSize(0.027);
  l2->SetHeader("MC");
  l2->AddEntry(gmc_mean1,"Photon 4 GeV, p_{T}^{leading #gamma} > 6 GeV","pe");
  l2->AddEntry(gmc_mean2,"MBD N&S>=1, p_{T}^{leading #gamma} > 1 GeV","pe");
  l2->AddEntry(gmc_mean3,"p_{T}/v_{z} weighted, p_{T}^{leading #gamma} > 1 GeV","pe");
  l2->Draw("same");
  c1->SaveAs("plot_pi0_peak_data_vs_mc.pdf");

  TCanvas* c2 = new TCanvas("c2","",700,700);
  c2->cd();
  c2->SetTicks(1,1);
  c2->SetRightMargin(0.10);
  c2->SetLeftMargin(0.13);
  c2->SetTopMargin(0.06);
  gdata_res1->Draw("AP");
  gdata_res2->Draw("P same");
  gmc_res1->Draw("P same");
  gmc_res2->Draw("P same");
  gdata_res1->GetXaxis()->SetTitle("p_{T}^{#pi^{0}} [GeV]");
  gdata_res1->GetYaxis()->SetTitle("m_{0} [GeV]");
  gdata_res1->GetYaxis()->CenterTitle();
  gdata_res1->GetXaxis()->CenterTitle();
  
  drawText("#bf{#it{sPHENIX}} Internal",sPHENIX_posx,sPHENIX_posy,1,22);
  drawText("pp #sqrt{s} = 200 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff,1,19);
  drawText("|v_{z}^{MBD}| < 30 cm",sPHENIX_posx,sPHENIX_posy-posy_diff*2,1,18);
  drawText("p_{T}^{#gamma} > 1 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff*3,1,18);
  drawText("#alpha < 0.6",sPHENIX_posx,sPHENIX_posy-posy_diff*4,1,18);
  drawText("#gamma prob. > 0.05",sPHENIX_posx,sPHENIX_posy-posy_diff*5,1,18);
  drawText("#pi^{0}",sPHENIX_posx,sPHENIX_posy-posy_diff*6,1,18);

  TLegend *l3 = new TLegend(0.18,0.71,0.37,0.88);
  SetLegendStyle(l3);
  l3->SetTextSize(0.027);
  l3->SetHeader("Data");
  l3->AddEntry(gdata_res1,"Photon 4 GeV, p_{T}^{leading #gamma} > 6 GeV","pe");
  l3->AddEntry(gdata_res2,"MBD N&S>=1, p_{T}^{leading #gamma} > 1 GeV","pe");
  l3->Draw("same");
  
  TLegend *l4 = new TLegend(0.18,0.52,0.37,0.69);
  SetLegendStyle(l4);
  l4->SetTextSize(0.027);
  l4->SetHeader("MC");
  l4->AddEntry(gmc_res1,"Photon 4 GeV, p_{T}^{leading #gamma} > 6 GeV","pe");
  l4->AddEntry(gmc_res2,"MBD N&S>=1, p_{T}^{leading #gamma} > 1 GeV","pe");
  l4->Draw("same");
  c2->SaveAs("plot_pi0_res_data_vs_mc.pdf");

}
