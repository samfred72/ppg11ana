#include <iostream>
#include <RooWorkspace.h>
#include <cmath>
#include "/sphenix/user/jpark4/sPHENIX_software/CaloTreeMaker/run/EmcalCluster/NeutralMesonAna/sPHENIX-PPG11Analysis/Analysis/Headers/ana.cxx"
#include "/sphenix/user/jpark4/sPHENIX_software/CaloTreeMaker/run/EmcalCluster/NeutralMesonAna/sPHENIX-PPG11Analysis/Analysis/Headers/TreeSetting.h"
#include "/sphenix/user/jpark4/Utility/Style_jaebeom.h"
#include "/sphenix/user/jpark4/Utility/commonUtility.h"

  
double GetGraphMaximum(TGraph* graph) {
  return graph ? TMath::MaxElement(graph->GetN(), graph->GetY()) : 0;
}

void CompCrossData(std::string particle="pi0",bool smear = 0)
{
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);
  ana anac;
  const char * smearstr = "";
  if (smear) smearstr="_smear";
  const int nbins = anac.nPtBins;
  const double* bins = anac.ptBins;
  double sPHENIX_posx = anac.sPHENIX_posx;
  double sPHENIX_posy = anac.sPHENIX_posy;
  double posy_diff = anac.posy_diff;
  
  double truthz_fraction = 2.7916; 
  TFile *f1 = new TFile(Form("/sphenix/user/samfred/run25/ppg11/histmaking/yields_without_prob_1GeV/data_%syields.root",particle.c_str()),"read");
  TH1D  *h1 = (TH1D*) f1->Get(Form("%syields",particle.c_str()));
  //h1->Scale(truthz_fraction);

  double lumidata = 0.03228*1e+09;
  double lumi = 200*1e6 / 42.0;
  double dy = 2.0;
  double br = (particle=="pi0") ? 0.98823 : 1;
  double norm = 2*M_PI;
  double mbdeff = 0.7; 

  TFile *f2 = new TFile(Form("/sphenix/user/samfred/run25/ppg11/CrossSectionClosure/eff_output_noreweight_pythiaMB_%s.root",particle.c_str()),"read");
  TH1D* h2 = (TH1D*) f2->Get("heff");
  
  std::string fitfilename = (particle == "pi0") ? "/sphenix/user/samfred/run25/ppg11/CrossSectionClosure/fit_pythiamb_cross_section_pi0.root" : "/sphenix/user/samfred/run25/ppg11/CrossSectionClosure/fit_pythiamb_cross_section_tsallis_eta.root";
  TFile *f_fit = new TFile(fitfilename.c_str(),"read");
  TF1* fitw = (TF1*) f_fit->Get("f_pTdNdpt");
  TF1* fit = (TF1*) f_fit->Get("fit_restricted");

  TFile *ft = new TFile(Form("/sphenix/user/samfred/run25/ppg11/CrossSectionClosure/hist_truth_noreweight_%s_MB.root",particle.c_str()),"read");
  TH1D* htruth = (TH1D*) ft->Get("heff_sp_pt_truth");
  htruth->Scale(truthz_fraction);

  TGraphErrors *gn = new TGraphErrors();
  TGraphErrors *grtruth = new TGraphErrors();
  TGraphErrors *gc = new TGraphErrors();

  int ip=0;
  for(int i=1; i <= h1->GetNbinsX(); ++i){
    float ptlow = h1->GetBinCenter(i) - h1->GetBinWidth(i)/2;
    float pthigh = h1->GetBinCenter(i) + h1->GetBinWidth(i)/2;
    float yield = h1->GetBinContent(i); 
    if(yield==0 || i<=2) continue;
    float yielderr = h1->GetBinError(i); 
    float ptpos = (ptlow + pthigh) / 2;
    float dpt = (pthigh-ptlow);
    
    double num = fitw->Integral(ptlow, pthigh);
    double den = fit->Integral(ptlow, pthigh);
    double meanpt = num/den; 
    
    double normalization = lumidata * dy * br * norm * dpt * meanpt * mbdeff;
    
    yield = yield/normalization;
    yielderr = yielderr/normalization;

    gn->SetPoint(ip,ptpos, yield);
    gn->SetPointError(ip, dpt/2, yielderr);
    
    double corr = h2->GetBinContent(i);
    double correrr = h2->GetBinError(i);
    gc->SetPoint(ip,ptpos,corr);
    gc->SetPointError(ip,dpt/2,correrr);

    ip++;
  }
  
  TGraphErrors * gr = getGraphRatio(gn,gc);

  ip=0;
  for(int i=1; i <= htruth->GetNbinsX(); ++i){
    float ptlow = htruth->GetBinCenter(i) - htruth->GetBinWidth(i)/2;
    float pthigh = htruth->GetBinCenter(i) + htruth->GetBinWidth(i)/2;
    float yield = htruth->GetBinContent(i); 
    if(yield==0 || i<=2) continue;
    float yielderr = htruth->GetBinError(i); 
    float ptpos = (ptlow + pthigh) / 2;
    float dpt = (pthigh-ptlow);
    
    double num = fitw->Integral(ptlow, pthigh);
    double den = fit->Integral(ptlow, pthigh);
    double meanpt = num/den; 
    
    double normalization = lumi * dy * br* norm * dpt * meanpt;
    
    yield = yield/normalization;
    yielderr = yielderr/normalization;

    grtruth->SetPoint(ip,ptpos, yield);
    grtruth->SetPointError(ip, dpt/2, yielderr);
    ip++;
  }
  grtruth->SetMarkerStyle(kFullDiamond);
  grtruth->SetMarkerSize(1.7);
  grtruth->SetMarkerColor(kGreen+2);
  grtruth->SetLineColor(kGreen+2);
  grtruth->SetFillColor(kGreen+2);
  gr->GetXaxis()->SetTitle("p_{T} [GeV]");
  gr->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}}[mb GeV^{-2}c^{3}]");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->CenterTitle();
  f_fit->Close();
  f2->Close();
  f1->Close();
  
  TGraphAsymmErrors *grphenix; 
  TFile *filefit = new TFile(Form("/sphenix/user/samfred/run25/ppg11/CrossSectionClosure/cross_section_%s_fit.root",particle.c_str()),"read");
  TF1 *fitc= (TF1*) filefit->Get("fit");
  fitc->SetNpx(1000);
  fitc->SetLineColor(kBlack);
  TGraphErrors *gratio_phenix = new TGraphErrors(); //ratio to phenix
  TGraphErrors *gratio_sphenix = new TGraphErrors();
  TDirectory *dir;
  if(particle=="pi0"){
    TFile *fphenix = new TFile("/sphenix/user/samfred/run25/ppg11/CrossSectionClosure/PHENIX_pi0_pp.root","read");
    dir = (TDirectory*) fphenix->Get("Figure 1d-2");
    grphenix = (TGraphAsymmErrors*) dir->Get("Graph1D_y2");
    TH1F* hist_e1= (TH1F*) dir->Get("Hist1D_y2_e1");
    TH1F* hist_e2= (TH1F*) dir->Get("Hist1D_y2_e2");

    for(int i=0; i<grphenix->GetN(); i++){
      float xpos = grphenix->GetPointX(i);
      float ypos = grphenix->GetPointY(i);

      float err1 = hist_e1->GetBinContent(i+1);
      float err2 = hist_e2->GetBinContent(i+1);

      float binlow = hist_e1->GetBinLowEdge(i+1);
      float binwidth = hist_e1->GetBinWidth(i+1);

      float xerrlow = xpos - binlow;
      float xerrhigh = binlow + binwidth - xpos;

      double toterr = sqrt(err1*err1 + err2*err2);
      grphenix->SetPoint(i,xpos,ypos);
      grphenix->SetPointError(i,xerrlow,xerrhigh,toterr,toterr);
    }

    for(int i=0; i<gr->GetN(); i++){
      float x = gr->GetPointX(i);
      float y = gr->GetPointY(i);
      float yerr = gr->GetErrorY(i);
      float yerrsys = gr->GetErrorY(i);
      float xerr = gr->GetErrorX(i);
      float xerrsys = gr->GetErrorX(i);
      float fit_integral = fitc->Integral(x-xerr, x+xerr) / (2*xerr); 
      //float fitval = fitc->Eval(x);//fit_integral;//fit->Eval(x);
      float fitval = fit_integral;//fit->Eval(x);
      float ratio = y/fitval;
      float ratioerr = yerr/fitval;
      gratio_sphenix->SetPoint(i,x,ratio);
      gratio_sphenix->SetPointError(i,xerr,ratioerr);
    }
    for(int i=0;i<grphenix->GetN(); i++){ 
      float x = grphenix->GetPointX(i);
      float y = grphenix->GetPointY(i);
      float yerr = grphenix->GetErrorY(i);
      float xerr = grphenix->GetErrorX(i);
      float fitval = fitc->Eval(x);
      float ratio = y/fitval;
      float ratioerr = yerr/fitval;
      gratio_phenix->SetPoint(i,x,ratio);
      gratio_phenix->SetPointError(i,xerr,ratioerr);
    }
  }
  else if(particle=="eta"){
    TFile *fphenix = new TFile("/sphenix/user/samfred/run25/ppg11/CrossSectionClosure/PHENIX_eta_pp.root","read");
    dir = (TDirectory*) fphenix->Get("Figure 2-10");
    grphenix = (TGraphAsymmErrors*) dir->Get("Graph1D_y1");
    TH1F* hist_eplus = (TH1F*) dir->Get("Hist1D_y1_e1plus");
    TH1F* hist_eminus = (TH1F*) dir->Get("Hist1D_y1_e1minus");
    TH1F* hist_e2= (TH1F*) dir->Get("Hist1D_y1_e2");
    TH1F* hist_e3= (TH1F*) dir->Get("Hist1D_y1_e3");
    TH1F* hist_e4= (TH1F*) dir->Get("Hist1D_y1_e4");

    for(int i=0; i<grphenix->GetN(); i++){
      float xpos = grphenix->GetPointX(i);
      float ypos = grphenix->GetPointY(i);
      float xerr = grphenix->GetErrorX(i);

      double errpl = hist_eplus->GetBinContent(i+1);
      double errmi = hist_eminus->GetBinContent(i+1);
      double err2 = hist_e2->GetBinContent(i+1);
      double err3 = hist_e3->GetBinContent(i+1);
      double err4 = hist_e4->GetBinContent(i+1);

      double toterr = sqrt(errpl*errpl + errmi*errmi + err2*err2 + err3*err3 + err4*err4);
      grphenix->SetPoint(i,xpos,ypos);
      grphenix->SetPointError(i,xerr,xerr,toterr,toterr);
    }
    int np1 = grphenix->GetN();

    TFile *fphenix2 = new TFile("/sphenix/user/samfred/run25/ppg11/CrossSectionClosure/PHENIX_eta_pp_v1.root","read");
    TDirectory *dir2 = (TDirectory*) fphenix2->Get("Figure 12.1");
    TGraphAsymmErrors *grphenix2 = (TGraphAsymmErrors*) dir2->Get("Graph1D_y1");
    TH1F* hist2_e1= (TH1F*) dir2->Get("Hist1D_y1_e1");
    TH1F* hist2_e2= (TH1F*) dir2->Get("Hist1D_y1_e2");
    TH1F* hist2_e3= (TH1F*) dir2->Get("Hist1D_y1_e3");
    TH1F* hist2_e4= (TH1F*) dir2->Get("Hist1D_y1_e4");

    for(int i=0; i<grphenix2->GetN(); i++){
      float xpos = grphenix2->GetPointX(i);
      float ypos = grphenix2->GetPointY(i);
      float xerr = grphenix2->GetErrorX(i);

      double err1 = hist2_e1->GetBinContent(i+1);
      double err2 = hist2_e2->GetBinContent(i+1);
      double err3 = hist2_e3->GetBinContent(i+1);
      double err4 = hist2_e4->GetBinContent(i+1);

      double toterr = sqrt(err1*err1 + err2*err2 + err3*err3 + err4*err4);
      grphenix->SetPoint(i + np1,xpos,ypos);
      grphenix->SetPointError(i + np1,xerr,xerr,toterr,toterr);
    }

    for(int i=0; i<gr->GetN(); i++){
      float x = gr->GetPointX(i);
      float y = gr->GetPointY(i);
      float yerr = gr->GetErrorY(i);
      float yerrsys = gr->GetErrorY(i);
      float xerr = gr->GetErrorX(i);
      float xerrsys = gr->GetErrorX(i);
      float fit_integral = fitc->Integral(x-xerr, x+xerr) / (2*xerr); 
      //float fitval = fitc->Eval(x);//fit_integral;//fit->Eval(x);
      float fitval = fit_integral;//fit->Eval(x);
      float ratio = y/fitval;
      float ratioerr = yerr/fitval;
      gratio_sphenix->SetPoint(i,x,ratio);
      gratio_sphenix->SetPointError(i,xerr,ratioerr);
    }
    for(int i=0;i<grphenix->GetN(); i++){ 
      float x = grphenix->GetPointX(i);
      float y = grphenix->GetPointY(i);
      float yerr = grphenix->GetErrorY(i);
      float xerr = grphenix->GetErrorX(i);
      float fitval = fitc->Eval(x);
      float ratio = y/fitval;
      float ratioerr = yerr/fitval;
      gratio_phenix->SetPoint(i,x,ratio);
      gratio_phenix->SetPointError(i,xerr,ratioerr);
    }

  }
  
  SetGraphStyle(grphenix,1,1);
  grphenix->SetMarkerColor(kGray+2);
  grphenix->SetLineColor(kBlack);
  SetGraphStyle(gr,0,0);
  grphenix->SetTitle(" ");
  grphenix->SetMinimum(5*1e-9);
  grphenix->SetMaximum(100);
  gr->SetMinimum(5*1e-9);
  gr->SetMaximum(100);
  grtruth->SetMinimum(5*1e-9);
  grtruth->SetMaximum(100);


  grphenix->GetXaxis()->SetLimits(0,10);
  gr->GetXaxis()->SetLimits(0,10);
  grtruth->GetXaxis()->SetLimits(0,10);
  grphenix->GetXaxis()->SetRangeUser(0,10);
  gr->GetXaxis()->SetRangeUser(0,10);
  grtruth->GetXaxis()->SetRangeUser(0,10);
  TCanvas* c1 = new TCanvas("c1","",700,700);
  TPad *pad1 = new TPad("pad1","",0,0.3,1,1);
  TPad *pad2 = new TPad("pad2","",0,0,1,0.3);
  c1->SetBottomMargin(10.22);
  //pad1->Draw();
  pad1->cd();
  pad1->SetTicks(1,1);
  pad1->SetRightMargin(0.09);
  pad1->SetLeftMargin(0.165);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.0);
  pad1->SetLogy();

  grphenix->GetYaxis()->SetTitleSize(0.041);
  grphenix->GetYaxis()->SetTitleOffset(1.9);
  grphenix->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} [mb GeV^{-2}c^{3}]");
  grphenix->GetYaxis()->CenterTitle();
  grphenix->GetXaxis()->CenterTitle();
  grphenix->GetXaxis()->SetTitleSize(0);
  grphenix->GetXaxis()->SetLabelSize(0);
  grphenix->GetYaxis()->SetLabelSize(0.05);
  grphenix->Draw("APE");
  gr->Draw("PE same");
  //grtruth->Draw("PE same");
  fitc->Draw("same");
  
  drawText("#bf{#it{sPHENIX}} Internal",sPHENIX_posx-0.03,sPHENIX_posy,1,22);
  drawText("#it{p+p} #kern[-0.5]{#sqrt{s}} = 200 GeV",sPHENIX_posx-0.03,sPHENIX_posy-posy_diff,1,19);

  TLegend *l = new TLegend(0.60,0.45,0.85,0.77);
  SetLegendStyle(l);
  l->SetTextSize(0.0385);
  if(particle=="pi0") l->SetHeader("#bf{#pi^{0}}");
  else l->SetHeader("#bf{#eta}");
  l->AddEntry(gr,"sPHENIX |y|<1","pe");
  l->AddEntry(grphenix,"PHENIX |y|<0.35","pe");
  //l->AddEntry(grtruth,"Pythia8 |y|<1","pe");
  l->AddEntry(fitc,"Fit to PHENIX","l");

  l->Draw("same");

  TLegend *l1 = new TLegend(0.51,0.52,0.65,0.62);
  SetLegendStyle(l1);
  l1->SetTextSize(0.0385);
  l1->AddEntry(fit,"Fit to PHENIX","l");
//  l1->Draw("same");
  c1->Update();

  SetGraphStyle(gratio_phenix,0,4);
  gratio_phenix->SetMarkerStyle(kFullSquare);
  gratio_phenix->SetMarkerColor(kGray+2);
  gratio_phenix->SetLineColor(kBlack);
  SetGraphStyle(gratio_sphenix,0,4);
  gratio_sphenix->SetMarkerStyle(kFullCircle);
  gratio_sphenix->SetMarkerColor(kPink-6);
  gratio_sphenix->SetLineColor(kPink-6);
//  pad2->Draw();
  pad2->cd();
  pad2->SetTicks(1,1);
  pad2->SetRightMargin(0.09);
  pad2->SetLeftMargin(0.165);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.32);
  pad2->SetLogy(0);
  gratio_phenix->GetXaxis()->SetLabelSize(0.11);
  gratio_phenix->GetYaxis()->SetLabelSize(0.11);
  gratio_phenix->GetYaxis()->SetNdivisions(505);
  gratio_phenix->GetXaxis()->SetTitleSize(0.12);
  gratio_phenix->GetYaxis()->SetTitleSize(0.11);
  gratio_phenix->GetYaxis()->SetTitle("Ratio-to-fit");
  if(particle=="pi0")  gratio_phenix->GetXaxis()->SetTitle("p_{T}^{#pi^{0}} [GeV/c]");
  else if(particle=="eta")  gratio_phenix->GetXaxis()->SetTitle("p_{T}^{#eta} [GeV/c]");
  gratio_phenix->GetXaxis()->CenterTitle();
  gratio_phenix->GetYaxis()->CenterTitle();
  gratio_phenix->GetYaxis()->SetTitleOffset(0.7);
  gratio_phenix->GetXaxis()->SetTitleOffset(1.1);
  gratio_phenix->GetYaxis()->SetLimits(0,3);
  gratio_phenix->GetYaxis()->SetRangeUser(0,2.99);
  gratio_phenix->GetXaxis()->SetLimits(0,10);
  gratio_phenix->GetXaxis()->SetRangeUser(0,10);
  gratio_sphenix->GetXaxis()->SetLimits(0,10);
  gratio_sphenix->GetXaxis()->SetRangeUser(0,10);

  gratio_phenix->Draw("AP");
  gratio_sphenix->Draw("P same");
  dashedLine(0,1,10,1,1,2);

  TLegend *l2 = new TLegend(0.56,0.78,0.88,0.96);
  SetLegendStyle(l2);
  l2->SetTextSize(0.075);
  l2->AddEntry(gratio_sphenix,"sPHENIX","pe");
  l2->AddEntry(gratio_phenix,"PHENIX","pe");
  l2->Draw("same");
  c1->cd();
  pad2->Draw();
  pad1->Draw();
  c1->Modified();
  c1->Update();
  c1->SaveAs(Form("/sphenix/user/samfred/run25/ppg11/CrossSectionClosure/plots/cross_data_%s.pdf",particle.c_str()));
}
