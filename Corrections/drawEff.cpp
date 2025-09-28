#include <iostream>
#include <cmath>
#include "../histmaking/ana.cxx"
#include "../histmaking/TreeSetting.h"
  
  
double GetGraphMaximum(TGraph* graph) {
  return graph ? TMath::MaxElement(graph->GetN(), graph->GetY()) : 0;
}

ana anac;
const int nbins = anac.nPtBins;
const double* bins = anac.ptBins;
double sPHENIX_posx = anac.sPHENIX_posx;
double sPHENIX_posy = anac.sPHENIX_posy;
double posy_diff = anac.posy_diff;


void drawEff(std::string particle="pi0",bool smear=0,const char * type = "MC")
{
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  std::string particlestr = particle=="pi0" ? "#pi^{0}" : "#eta";
  const char * smearstr = "";
  if (smear) smearstr = "_smear";
  TFile *f1 = new TFile(Form("/sphenix/user/samfred/run25/ppg11/CrossSectionClosure/hist_truth_noreweight_%s_MB.root",particle.c_str()),"read");
  TFile *f2 = new TFile(Form("/sphenix/user/samfred/run25/ppg11/histmaking/yields_without_prob_1GeV%s/%s_%syields.root",smearstr,type,particle.c_str()),"read");

  TH1D  *eff_den = (TH1D*) f1->Get("heff_sp_pt_den");
  TH1D  *eff_num = (TH1D*) f2->Get(Form("%syields",particle.c_str()));
  if (!eff_den || !eff_num) {
    std::cerr << "Error: One of the TEfficiency objects is null!" << std::endl;
    return;
  }
  for (int i = 1; i < 8; i++) {
    cout << eff_num->GetBinContent(i) << " " << eff_den->GetBinContent(i) << " " << eff_num->GetBinContent(i)/eff_den->GetBinContent(i) << endl;
  }

  int xrange = 20;

  TH1D* heff = (TH1D*) eff_num -> Clone("heff");
  cout << heff->GetNbinsX() << " " << eff_den->GetNbinsX() << endl;
  heff->Divide(heff,eff_den,1,1,"B");

  heff->GetXaxis()->SetRangeUser(0,xrange);
  heff->GetXaxis()->SetRangeUser(0,xrange);

  heff->GetYaxis()->SetLimits(0,1);
  heff->GetYaxis()->SetRangeUser(0,1);
  heff->GetXaxis()->SetTitle("p_{T}^{reco} [GeV/c]");
  heff->GetYaxis()->SetTitle("Acceptance x Efficiency");
  heff->GetXaxis()->CenterTitle();
  heff->GetYaxis()->CenterTitle();
  
  heff->SetLineColor(kRed+1);
  heff->SetMarkerColor(kRed+1);

  heff->SetLineWidth(5);
  
  TCanvas *c1 = new TCanvas("c1","",700,700);
  c1->SetTicks(1,1);
  c1->SetRightMargin(0.10);
  c1->SetLeftMargin(0.13);
  c1->SetTopMargin(0.06);
  c1->cd();
  heff->Draw("PE");
  TLegend *l = new TLegend(0.18,0.88,0.4,0.59);
  SetLegendStyle(l);
  l->SetHeader(Form("#bf{%s} #rightarrow #gamma+#gamma",particlestr.c_str()));
  l->AddEntry(heff,"Pythia8 MB","le");
  l->Draw("same");

  drawText("#bf{#it{sPHENIX}} Internal",sPHENIX_posx,sPHENIX_posy,1,22);
  drawText("pp #sqrt{s} = 200 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff,1,19);
  drawText("|v_{z}^{truth}| < 30 cm",sPHENIX_posx,sPHENIX_posy-posy_diff*2,1,18);
  drawText("p_{T}^{#gamma} > 1 GeV",sPHENIX_posx,sPHENIX_posy-posy_diff*3,1,18);
  drawText("|#eta^{truth}| < 1",sPHENIX_posx,sPHENIX_posy-posy_diff*4,1,18);
  drawText("#alpha < 0.6",sPHENIX_posx,sPHENIX_posy-posy_diff*5,1,18);
  //drawText("#gamma prob. > 0.05",sPHENIX_posx,sPHENIX_posy-posy_diff*6,1,18);

  c1->SaveAs(Form("/sphenix/user/samfred/run25/ppg11/CrossSectionClosure/eff_noweight_%s_pythia.pdf",particle.c_str()));


  TFile *wf = new TFile(Form("/sphenix/user/samfred/run25/ppg11/CrossSectionClosure/eff_output_noreweight_pythiaMB_%s.root",particle.c_str()),"recreate");
  wf->cd();
  heff->Write();
}
