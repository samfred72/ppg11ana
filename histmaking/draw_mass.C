#include "../Headers/Style_jaebeom.h"
#include "../Headers/sPhenixStyle.h"

void draw_mass(const char * type = "Jet10") {
  TFile * f = TFile::Open(Form("mass_%s.root",type),"READ");
  TH1D * masshists[25];
  TCanvas * c = new TCanvas("c","",500,500);
  c->SaveAs(Form("mass_%s_hists.pdf[",type));
  gStyle->SetOptStat(0);
  for (int i = 0; i < 25; i++) {
    masshists[i] = (TH1D*)f->Get(Form("hmass_pt%i",i));
    masshists[i]->SetLineColor(1);
    masshists[i]->SetLineWidth(2);
    masshists[i]->Draw();
    drawText("#bf{#it{sPHENIX}} Internal",0.35,0.81,1,22);
    drawText("MC Pythia run 21",0.35,0.74,1,22);
    drawText(Form("%i < pT < %i",i,i+1),0.5,0.68,1,16);
    drawText("|vz| < 30 cm",0.5,0.63,1,16);
    drawText("|#eta| < 1.1",0.5,0.58,1,16);
    drawText("p_{T,#gamma} > 0.5 GeV",0.5,0.53,1,16);
    drawText("prob_{#gamma} > 0.05",0.5,0.48,1,16);
    drawText("#alpha < 0.6",0.5,0.43,1,16);

    c->SaveAs(Form("mass_%s_hists.pdf",type));
  }
  c->SaveAs(Form("mass_%s_hists.pdf]",type));
  return;
}

