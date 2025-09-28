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

void drawdata(const char * particle="pi0")
{
  const char * formattedparticle = "#pi^{0}";
  if (strcmp(particle,"eta") == 0) formattedparticle = "#eta";
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);
  ana anac;
  const int nbins = anac.nPtBins;
  const double* bins = anac.ptBins;
  double sPHENIX_posx = anac.sPHENIX_posx;
  double sPHENIX_posy = anac.sPHENIX_posy;
  double posy_diff = anac.posy_diff;

  TGraphErrors* gdata_mean = new TGraphErrors();
  TGraphErrors* gdata_res = new TGraphErrors();
  TGraphErrors* gdata_ratio = new TGraphErrors();
  TGraphErrors* gmc_mean = new TGraphErrors();
  TGraphErrors* gmc_res = new TGraphErrors();
  TGraphErrors* gmc_ratio = new TGraphErrors();
  TGraphErrors* gmcs_mean = new TGraphErrors();
  TGraphErrors* gmcs_res = new TGraphErrors();
  TGraphErrors* gmcs_ratio = new TGraphErrors();
  TGraphErrors* gmcs2_mean = new TGraphErrors();
  TGraphErrors* gmcs2_res = new TGraphErrors();
  TGraphErrors* gmcs2_ratio = new TGraphErrors();

  int pointcount=0;
  int datapointcount = 0;
  int mcpointcount = 0;
  int smearpointcount = 0;
  for(int i = 2; i < 8; ++i) {
    const char * filename;
    RooWorkspace *w;
    RooFitResult *fitRes;
    RooRealVar *mean;
    RooRealVar *sigma;
    RooFormulaVar *ratio;
    TFile * f;
    filename = Form("/sphenix/user/samfred/run25/ppg11/histmaking/sufficient_without_prob_1GeV/workspace_fits_data_%s_pt%i_mbd.root",particle,i);
    if (!gSystem->AccessPathName(filename)) {
      f = TFile::Open(filename,"READ");
      w = (RooWorkspace*) f->Get("workspace");
      fitRes = (RooFitResult*)w->obj("fitresult_model_binnedData");

      //DATA
      mean = w->var("mean"); 
      sigma = w->var("sigma1"); 
      ratio = new RooFormulaVar("ratio", "sigma1 / mean", RooArgList(*sigma, *mean));

      gdata_mean->SetPoint(datapointcount, (i+0.5), mean->getVal());
      gdata_mean->SetPointError(datapointcount, 0.5, mean->getError());

      gdata_res->SetPoint(datapointcount, (i+0.5), sigma->getVal());
      gdata_res->SetPointError(datapointcount, 0.5, sigma->getError());

      gdata_ratio->SetPoint(datapointcount, (i+0.5), ratio->getVal());
      gdata_ratio->SetPointError(datapointcount, 0.5, ratio->getPropagatedError(*fitRes));

      f->Close();
      delete f;
      datapointcount++;
    }
    filename = Form("/sphenix/user/samfred/run25/ppg11/histmaking/sufficient_without_prob_1GeV/workspace_fits_MC_%s_pt%i_MB.root",particle,i);
    if (!gSystem->AccessPathName(filename)) {
      f = TFile::Open(filename,"READ");
      w = (RooWorkspace*) f->Get("workspace");
      fitRes = (RooFitResult*)w->obj("fitresult_model_binnedData");

      //MC
      mean = w->var("mean"); 
      sigma = w->var("sigma1"); 
      ratio = new RooFormulaVar("ratio", "sigma1 / mean", RooArgList(*sigma, *mean));

      gmc_mean->SetPoint(mcpointcount, (i+0.5), mean->getVal());
      gmc_mean->SetPointError(mcpointcount, 0.5, mean->getError());

      gmc_res->SetPoint(mcpointcount, (i+0.5), sigma->getVal());
      gmc_res->SetPointError(mcpointcount, 0.5, sigma->getError());

      gmc_ratio->SetPoint(mcpointcount, (i+0.5), ratio->getVal());
      gmc_ratio->SetPointError(mcpointcount, 0.5, ratio->getPropagatedError(*fitRes));

      f->Close();
      delete f;
      mcpointcount++;
    }
    filename = Form("/sphenix/user/samfred/run25/ppg11/histmaking/sufficient_without_prob_1GeV_smear/workspace_fits_MC_%s_pt%i_MB.root",particle,i);
    if (!gSystem->AccessPathName(filename)) {
      f = TFile::Open(filename,"READ");
      w = (RooWorkspace*) f->Get("workspace");
      fitRes = (RooFitResult*)w->obj("fitresult_model_binnedData");

      //MC
      mean = w->var("mean"); 
      sigma = w->var("sigma1"); 
      ratio = new RooFormulaVar("ratio", "sigma1 / mean", RooArgList(*sigma, *mean));

      gmcs_mean->SetPoint(smearpointcount, (i+0.5), mean->getVal());
      gmcs_mean->SetPointError(smearpointcount, 0.5, mean->getError());

      gmcs_res->SetPoint(smearpointcount, (i+0.5), sigma->getVal());
      gmcs_res->SetPointError(smearpointcount, 0.5, sigma->getError());

      gmcs_ratio->SetPoint(smearpointcount, (i+0.5), ratio->getVal());
      gmcs_ratio->SetPointError(smearpointcount, 0.5, ratio->getPropagatedError(*fitRes));


      f->Close();
      delete f;
      smearpointcount++;
    }
    /*
       f = new TFile(Form("/sphenix/user/samfred/run25/ppg11/histmaking/sufficient/workspace_fits_MC_%s_pt%i.root",particle,i),"read");

       w = (RooWorkspace*) f->Get("workspace");
       fitRes = (RooFitResult*)w->obj("fitresult_model_binnedData");

    //MC
    mean = w->var("mean"); 
    sigma = w->var("sigma1"); 
    ratio = new RooFormulaVar("ratio", "sigma1 / mean", RooArgList(*sigma, *mean));

    gmcs2_mean->SetPoint(pointcount, (i+0.5), mean->getVal());
    gmcs2_mean->SetPointError(pointcount, 0.5, mean->getError());

    gmcs2_res->SetPoint(pointcount, (i+0.5), sigma->getVal());
    gmcs2_res->SetPointError(pointcount, 0.5, sigma->getError());

    gmcs2_ratio->SetPoint(pointcount, (i+0.5), ratio->getVal());
    gmcs2_ratio->SetPointError(pointcount, 0.5, ratio->getPropagatedError(*fitRes));


    f->Close();
    */  
    //pointcount++;
  }
  TGraphErrors * gratio = getGraphRatio(gdata_mean,gmc_mean);
  
  gdata_mean->GetXaxis()->SetLimits(0,10);
  gdata_mean->GetXaxis()->SetRangeUser(0,10);
  gdata_mean->GetXaxis()->SetLimits(0,10);
  gdata_mean->GetXaxis()->SetRangeUser(0,10);
  gdata_res->GetXaxis()->SetLimits(0,10);
  gdata_res->GetXaxis()->SetRangeUser(0,10);
  gdata_res->GetXaxis()->SetLimits(0,10);
  gdata_res->GetXaxis()->SetRangeUser(0,10);
  gdata_ratio->GetXaxis()->SetLimits(0,10);
  gdata_ratio->GetXaxis()->SetRangeUser(0,10);
  gdata_ratio->GetXaxis()->SetLimits(0,10);
  gdata_ratio->GetXaxis()->SetRangeUser(0,10);
  gratio->GetXaxis()->SetLimits(0,10);
  gratio->GetXaxis()->SetRangeUser(0,10);

  SetGraphStyle(gdata_mean,0,0);
  SetGraphStyle(gdata_res,0,0);
  SetGraphStyle(gdata_ratio,0,0);
  SetGraphStyle(gratio,0,0);

  //gdata_mean->SetMarkerStyle(kOpenCircle);
  //gdata_res->SetMarkerStyle(kOpenCircle);
  //gdata_ratio->SetMarkerStyle(kOpenCircle);
 
  float limitmeanlow =   (strcmp("pi0",particle) == 0) ? 0.05 : 0.30;
  float limitmeanhigh =  (strcmp("pi0",particle) == 0) ? 0.30 : 0.70;
  float rangemeanlow =   (strcmp("pi0",particle) == 0) ? 0.11 : 0.50;
  float rangemeanhigh =  (strcmp("pi0",particle) == 0) ? 0.20 : 0.75;
  float limitreslow =    (strcmp("pi0",particle) == 0) ? 0.00 : 0.00;
  float limitreshigh =   (strcmp("pi0",particle) == 0) ? 0.05 : 0.20;
  float rangereslow =    (strcmp("pi0",particle) == 0) ? 0.00 : 0.00;
  float rangereshigh =   (strcmp("pi0",particle) == 0) ? 0.05 : 0.20;
  float limitratiolow =  (strcmp("pi0",particle) == 0) ? 0.00 : 0.00;
  float limitratiohigh = (strcmp("pi0",particle) == 0) ? 0.50 : 0.50;
  float rangeratiolow =  (strcmp("pi0",particle) == 0) ? 0.08 : 0.02;
  float rangeratiohigh = (strcmp("pi0",particle) == 0) ? 0.30 : 0.22;
  float limitratlow =    (strcmp("pi0",particle) == 0) ? 0.80 : 0.80;
  float limitrathigh =   (strcmp("pi0",particle) == 0) ? 1.50 : 1.50;
  float rangeratlow =    (strcmp("pi0",particle) == 0) ? 0.80 : 0.80;
  float rangerathigh =   (strcmp("pi0",particle) == 0) ? 1.50 : 1.50;
  gdata_mean->GetYaxis()->SetLimits(    limitmeanlow,limitmeanhigh);
  gdata_mean->GetYaxis()->SetRangeUser( rangemeanlow,rangemeanhigh);
  gdata_res->GetYaxis()->SetLimits(     limitreslow,limitreshigh);
  gdata_res->GetYaxis()->SetRangeUser(  rangereslow,rangereshigh);
  gdata_ratio->GetYaxis()->SetLimits(   limitratiolow,limitratiohigh);
  gdata_ratio->GetYaxis()->SetRangeUser(rangeratiolow,rangeratiohigh);
  gratio->GetYaxis()->SetLimits(limitratlow,limitrathigh);
  gratio->GetYaxis()->SetRangeUser(rangeratlow,rangerathigh);

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
  
  gmc_mean->GetYaxis()->SetLimits(    limitmeanlow,limitmeanhigh);
  gmc_mean->GetYaxis()->SetRangeUser( rangemeanlow,rangemeanhigh);
  gmc_res->GetYaxis()->SetLimits(     limitreslow,limitreshigh);
  gmc_res->GetYaxis()->SetRangeUser(  rangereslow,rangereshigh);
  gmc_ratio->GetYaxis()->SetLimits(   limitratiolow,limitratiohigh);
  gmc_ratio->GetYaxis()->SetRangeUser(rangeratiolow,rangeratiohigh);

  gmcs_mean->GetXaxis()->SetLimits(0,10);
  gmcs_mean->GetXaxis()->SetRangeUser(0,10);
  gmcs_mean->GetXaxis()->SetLimits(0,10);
  gmcs_mean->GetXaxis()->SetRangeUser(0,10);
  gmcs_res->GetXaxis()->SetLimits(0,10);
  gmcs_res->GetXaxis()->SetRangeUser(0,10);
  gmcs_res->GetXaxis()->SetLimits(0,10);
  gmcs_res->GetXaxis()->SetRangeUser(0,10);
  gmcs_ratio->GetXaxis()->SetLimits(0,10);
  gmcs_ratio->GetXaxis()->SetRangeUser(0,10);
  gmcs_ratio->GetXaxis()->SetLimits(0,10);
  gmcs_ratio->GetXaxis()->SetRangeUser(0,10);

  SetGraphStyle(gmcs_mean,1,0);
  SetGraphStyle(gmcs_res,1,0);
  SetGraphStyle(gmcs_ratio,1,0);

  gmcs_mean->SetMarkerStyle(kOpenCircle);
  gmcs_res->SetMarkerStyle(kOpenCircle);
  gmcs_ratio->SetMarkerStyle(kOpenCircle);
  
  gmcs_mean->GetYaxis()->SetLimits(    limitmeanlow,limitmeanhigh);
  gmcs_mean->GetYaxis()->SetRangeUser( rangemeanlow,rangemeanhigh);
  gmcs_res->GetYaxis()->SetLimits(     limitreslow,limitreshigh);
  gmcs_res->GetYaxis()->SetRangeUser(  rangereslow,rangereshigh);
  gmcs_ratio->GetYaxis()->SetLimits(   limitratiolow,limitratiohigh);
  gmcs_ratio->GetYaxis()->SetRangeUser(rangeratiolow,rangeratiohigh);

  gmcs2_mean->GetXaxis()->SetLimits(0,10);
  gmcs2_mean->GetXaxis()->SetRangeUser(0,10);
  gmcs2_mean->GetXaxis()->SetLimits(0,10);
  gmcs2_mean->GetXaxis()->SetRangeUser(0,10);
  gmcs2_res->GetXaxis()->SetLimits(0,10);
  gmcs2_res->GetXaxis()->SetRangeUser(0,10);
  gmcs2_res->GetXaxis()->SetLimits(0,10);
  gmcs2_res->GetXaxis()->SetRangeUser(0,10);
  gmcs2_ratio->GetXaxis()->SetLimits(0,10);
  gmcs2_ratio->GetXaxis()->SetRangeUser(0,10);
  gmcs2_ratio->GetXaxis()->SetLimits(0,10);
  gmcs2_ratio->GetXaxis()->SetRangeUser(0,10);

  SetGraphStyle(gmcs2_mean,2,0);
  SetGraphStyle(gmcs2_res,2,0);
  SetGraphStyle(gmcs2_ratio,2,0);

  gmcs2_mean->SetMarkerStyle(kOpenCircle);
  gmcs2_res->SetMarkerStyle(kOpenCircle);
  gmcs2_ratio->SetMarkerStyle(kOpenCircle);
  
  gmcs2_mean->GetYaxis()->SetLimits(    limitmeanlow,limitmeanhigh);
  gmcs2_mean->GetYaxis()->SetRangeUser( rangemeanlow,rangemeanhigh);
  gmcs2_res->GetYaxis()->SetLimits(     limitreslow,limitreshigh);
  gmcs2_res->GetYaxis()->SetRangeUser(  rangereslow,rangereshigh);
  gmcs2_ratio->GetYaxis()->SetLimits(   limitratiolow,limitratiohigh);
  gmcs2_ratio->GetYaxis()->SetRangeUser(rangeratiolow,rangeratiohigh);
  
  // mean
  TCanvas* c1 = new TCanvas("c1","",700,700);
  c1->cd();
  c1->SetTicks(1,1);
  c1->SetRightMargin(0.10);
  c1->SetLeftMargin(0.13);
  c1->SetTopMargin(0.06);
  gmc_mean->Draw("AP");
  gmc_mean->GetXaxis()->SetTitle(Form("p_{T}^{%s} [GeV]",formattedparticle));
  gmc_mean->GetYaxis()->SetTitle("m_{0} [GeV]");
  gmc_mean->GetYaxis()->CenterTitle();
  gmc_mean->GetXaxis()->CenterTitle();
  gdata_mean->Draw("same P");
  gmcs_mean->Draw("same P");
  //gmcs2_mean->Draw("same P");
  
  drawText("#bf{#it{sPHENIX}} Internal",sPHENIX_posx,sPHENIX_posy,1,22);
  drawText("pp #sqrt{s} = 200 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff,1,19);
  drawText("|v_{z}^{MBD}| < 30 cm",sPHENIX_posx,sPHENIX_posy-posy_diff*2,1,18);
  drawText("p_{T}^{#gamma} > 1 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff*3,1,18);
  drawText("#alpha < 0.6",sPHENIX_posx,sPHENIX_posy-posy_diff*4,1,18);
  //drawText("#gamma prob. > 0.05",sPHENIX_posx,sPHENIX_posy-posy_diff*5,1,18);
  drawText(formattedparticle,sPHENIX_posx,sPHENIX_posy-posy_diff*6,1,18);

  TLegend *l1 = new TLegend(0.18,0.71,0.37,0.88);
  SetLegendStyle(l1);
  l1->SetTextSize(0.027);
  l1->SetHeader("Run 47289-53879");
  l1->AddEntry(gdata_mean,"MB scaled, p_{T}^{leading #gamma} > 1 GeV","pe");
  l1->AddEntry(gmc_mean,"Pythia MB, p_{T}^{leading #gamma} > 1 GeV","pe");
  //l1->AddEntry(gmcs_mean,"Pythia MB smeared 0.001, p_{T}^{leading #gamma} > 1 GeV","pe");
  l1->AddEntry(gmcs_mean,"Pythia MB smeared, p_{T}^{leading #gamma} > 1 GeV","pe");
  l1->Draw("same");
  
  //c1->SaveAs(Form("plot_%s_peak_data.pdf",particle));

  // width
  TCanvas* c2 = new TCanvas("c2","",700,700);
  c2->cd();
  c2->SetTicks(1,1);
  c2->SetRightMargin(0.10);
  c2->SetLeftMargin(0.13);
  c2->SetTopMargin(0.06);
  gmc_res->Draw("AP");
  gmc_res->GetXaxis()->SetTitle(Form("p_{T}^{%s} [GeV]",formattedparticle));
  gmc_res->GetYaxis()->SetTitle("m_{0} width [GeV]");
  gmc_res->GetYaxis()->CenterTitle();
  gmc_res->GetXaxis()->CenterTitle();
  gdata_res->Draw("same P");
  gmcs_res->Draw("same P");
  gmcs2_res->Draw("same P");
  
  drawText("#bf{#it{sPHENIX}} Internal",sPHENIX_posx,sPHENIX_posy,1,22);
  drawText("pp #sqrt{s} = 200 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff,1,19);
  drawText("|v_{z}^{MBD}| < 30 cm",sPHENIX_posx,sPHENIX_posy-posy_diff*2,1,18);
  drawText("p_{T}^{#gamma} > 1 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff*3,1,18);
  drawText("#alpha < 0.6",sPHENIX_posx,sPHENIX_posy-posy_diff*4,1,18);
  drawText(formattedparticle,sPHENIX_posx,sPHENIX_posy-posy_diff*5,1,18);
  //drawText("#gamma prob. > 0.05",sPHENIX_posx,sPHENIX_posy-posy_diff*6,1,18);
  
  TLegend *l2 = new TLegend(0.18,0.71,0.37,0.88);
  SetLegendStyle(l2);
  l2->SetTextSize(0.027);
  l2->SetHeader("Run 47289-53879");
  l2->AddEntry(gdata_res,"MB scaled, p_{T}^{leading #gamma} > 1 GeV","pe");
  l2->AddEntry(gmc_res,"Pythia MB, p_{T}^{leading #gamma} > 1 GeV","pe");
  //l2->AddEntry(gmcs_res,"Pythia MB smeared 0.001, p_{T}^{leading #gamma} > 1 GeV","pe");
  l2->AddEntry(gmcs_res,"Pythia MB smeared, p_{T}^{leading #gamma} > 1 GeV","pe");
  l2->Draw("same");
  

  // ratio
  TCanvas* c3 = new TCanvas("c3","",700,700);
  c3->cd();
  c3->SetTicks(1,1);
  c3->SetRightMargin(0.10);
  c3->SetLeftMargin(0.13);
  c3->SetTopMargin(0.06);
  gmc_ratio->Draw("AP");
  gmc_ratio->GetXaxis()->SetTitle(Form("p_{T}^{%s} [GeV]",formattedparticle));
  gmc_ratio->GetYaxis()->SetTitle("width/m_{0}");
  gmc_ratio->GetYaxis()->CenterTitle();
  gmc_ratio->GetXaxis()->CenterTitle();
  gdata_ratio->Draw("same P");
  gmcs_ratio->Draw("same P");
  //gmcs2_ratio->Draw("same P");
  
  drawText("#bf{#it{sPHENIX}} Internal",sPHENIX_posx,sPHENIX_posy,1,22);
  drawText("pp #sqrt{s} = 200 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff,1,19);
  drawText("|v_{z}^{MBD}| < 30 cm",sPHENIX_posx,sPHENIX_posy-posy_diff*2,1,18);
  drawText("p_{T}^{#gamma} > 1 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff*3,1,18);
  drawText("#alpha < 0.6",sPHENIX_posx,sPHENIX_posy-posy_diff*4,1,18);
  drawText(formattedparticle,sPHENIX_posx,sPHENIX_posy-posy_diff*5,1,18);
  
  TLegend *l3 = new TLegend(0.18,0.71,0.37,0.88);
  SetLegendStyle(l3);
  l3->SetTextSize(0.027);
  l3->SetHeader("Run 47289-53879");
  l3->AddEntry(gdata_res,"MB scaled, p_{T}^{leading #gamma} > 1 GeV","pe");
  l3->AddEntry(gmc_res,"Pythia MB, p_{T}^{leading #gamma} > 1 GeV","pe");
  //l3->AddEntry(gmcs_res,"Pythia MB smeared 0.001, p_{T}^{leading #gamma} > 1 GeV","pe");
  l3->AddEntry(gmcs_res,"Pythia MB smeared, p_{T}^{leading #gamma} > 1 GeV","pe");
  l3->Draw("same");
/*
  // data-mc ratio
  TCanvas* c4 = new TCanvas("c4","",700,700);
  c4->cd();
  c4->SetTicks(1,1);
  c4->SetRightMargin(0.10);
  c4->SetLeftMargin(0.13);
  c4->SetTopMargin(0.06);
  gratio->Draw("AP");
  gratio->GetXaxis()->SetTitle(Form("p_{T}^{%s} [GeV]",formattedparticle));
  gratio->GetYaxis()->SetTitle("m_{0}^{data}/m_{0}^{MC}");
  gratio->GetYaxis()->CenterTitle();
  gratio->GetXaxis()->CenterTitle();
  
  drawText("#bf{#it{sPHENIX}} Internal",sPHENIX_posx,sPHENIX_posy,1,22);
  drawText("pp #sqrt{s} = 200 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff,1,19);
  drawText("|v_{z}^{MBD}| < 30 cm",sPHENIX_posx,sPHENIX_posy-posy_diff*2,1,18);
  drawText("p_{T}^{#gamma} > 1 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff*3,1,18);
  drawText("#alpha < 0.6",sPHENIX_posx,sPHENIX_posy-posy_diff*4,1,18);
  drawText(formattedparticle,sPHENIX_posx,sPHENIX_posy-posy_diff*5,1,18);
  
  TLegend *l4 = new TLegend(0.18,0.71,0.37,0.88);
  SetLegendStyle(l4);
  l4->SetTextSize(0.027);
  l4->SetHeader("Run 47289-53879");
  l4->AddEntry(gratio,"DATA/MC, p_{T}^{leading #gamma} > 1 GeV","pe");
  l4->Draw("same");

  TLine * l = new TLine(0,1,10,1);
  l->SetLineStyle(kDashed);
  l->Draw("same");
*/
}
