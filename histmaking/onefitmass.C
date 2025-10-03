#include "commonUtility.h"
#include "Style_jaebeom.h"
#include "onefitmass.h"
#include <map>

using namespace RooFit;

const char * workspacename;
const char * workspacesufficient;
const char * pdffilename;
const char * pdfsufficient;

void move() {
  // Define the source and destination paths
  std::filesystem::path source_path = workspacename;
  std::filesystem::path destination_path = workspacesufficient;

  try {
    // Attempt to move the file
    std::filesystem::rename(source_path, destination_path);
    std::cout << "File moved successfully from " << source_path << " to " << destination_path << std::endl;
  } catch (const std::filesystem::filesystem_error& e) {
    // Handle potential errors during the move operation
    std::cerr << "Error moving file: " << e.what() << std::endl;
  }
  
  source_path = pdffilename;
  destination_path = pdfsufficient;

  try {
    // Attempt to move the file
    std::filesystem::rename(source_path, destination_path);
    std::cout << "File moved successfully from " << source_path << " to " << destination_path << std::endl;
  } catch (const std::filesystem::filesystem_error& e) {
    // Handle potential errors during the move operation
    std::cerr << "Error moving file: " << e.what() << std::endl;
  }
  return;
}

TH1D * shorthist(TH1D * hist, float low, float high) {
  int teststart = hist->FindBin(low);
  int zerostart = 0;
  for (int i = 0; i < hist->GetNbinsX(); i++) {
    if (hist->GetBinContent(i) == 0) zerostart++;
  }
  int binstart = teststart;// max(teststart,zerostart);
  int binend = hist->FindBin(high);
  cout << "binstart: " << binstart << " binend: " << binend << endl;
  if (binend - binstart <= 0) return hist;
  float binwidth = hist->GetBinWidth(binstart);
  TH1D * returnhist = new TH1D(Form("%s_short",hist->GetName()),"",binend-binstart,hist->GetXaxis()->GetBinLowEdge(binstart),hist->GetXaxis()->GetBinUpEdge(binend));
  for (int i = binstart; i <= binend; i++) {
    returnhist->SetBinContent(i-binstart,hist->GetBinContent(i));
    returnhist->SetBinError(i-binstart,hist->GetBinError(i));
  }
  return returnhist;
}


void onefitmass(const char * particle = "pi0", const char * sample = "MC", const char * trigger = "MB", int pt = 1) {
  set_maps();

  map<string,int> textmap = {
                             {"pi0",0},{"eta",1},
                             {"data",0},{"MC",1},
                             {"MB",0},{"photon3",1},{"photon4",2},{"photon5",3},
                             {"Jet10",1},{"Jet20",2},{"Jet30",3}
                            };

  float plotlow =      map_plotlow     [textmap[particle]][textmap[sample]][textmap[trigger]][pt];
  float plothigh =     map_plothigh    [textmap[particle]][textmap[sample]][textmap[trigger]][pt];
  float initialmean =  map_initialmean [textmap[particle]][textmap[sample]][textmap[trigger]][pt];
  float initialsigma = map_initialsigma[textmap[particle]][textmap[sample]][textmap[trigger]][pt];
  float Alow =         map_Alow        [textmap[particle]][textmap[sample]][textmap[trigger]][pt];
  float Ahigh =        map_Ahigh       [textmap[particle]][textmap[sample]][textmap[trigger]][pt];
  float Astart =       map_Astart      [textmap[particle]][textmap[sample]][textmap[trigger]][pt];
  float Blow =         map_Blow        [textmap[particle]][textmap[sample]][textmap[trigger]][pt];
  float Bhigh =        map_Bhigh       [textmap[particle]][textmap[sample]][textmap[trigger]][pt];
  float Bstart =       map_Bstart      [textmap[particle]][textmap[sample]][textmap[trigger]][pt];

  TFile * wf = new TFile(Form("workspaces/workspace_fits_%s_%s_pt%i_%s.root",sample,particle,pt,trigger),"RECREATE"); 
  TFile * inf;
  if (strcmp(sample,"data") == 0) inf = TFile::Open("hists/mass_data.root","READ");
  else if (strcmp(sample,"MC") == 0) inf = TFile::Open(Form("hists/mass_%s_%s.root",sample,trigger),"READ");
  TH1D * h;
  if (strcmp(sample,"data") == 0) h = shorthist((TH1D*)inf->Get(Form("hmass_%s_pt%i_%s",particle,pt,trigger)),plotlow,plothigh);
  else if (strcmp(sample,"MC") == 0) h = shorthist((TH1D*)inf->Get(Form("hmass_%s_pt%i",particle,pt)),plotlow,plothigh);

  TCanvas * c = new TCanvas("c","",600,800);
  gStyle->SetOptStat(0);

  int binlow = h->FindBin(plotlow);
  int binhigh = h->FindBin(plothigh);
  int nbins = binhigh - binlow;
  cout << "Low bin: " << binlow << " High bin: " << binhigh << " nbins: " << nbins << endl;
  //RooMsgService::instance().setSilentMode(true);
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->factory("mass[0.0, 1.0]");
  ws->var("mass")->setRange(plotlow,plothigh);
  ws->var("mass")->Print();
  ws->var("mass")->setBins(nbins);

  float integral = h->Integral(1,h->GetNbinsX());
  if (integral < 10) {
    cout << "Histogram too empty! Exiting." << endl;
  }


  if (strcmp(particle,"eta") == 0 && strcmp(sample,"MC") == 0) {
    h->Rebin(2);
    //if (pt == 2) h->Rebin(2);
    if (pt > 3) h->Rebin(2);
    if (pt > 4) h->Rebin(2);
  }
  if (strcmp(particle,"pi0") == 0 && strcmp(sample,"MC") == 0) {
    if (pt >= 5) h->Rebin(2);
  }
  if (strcmp(particle,"pi0") == 0 && strcmp(sample,"data") == 0) {
    if (pt >= 7) h->Rebin(2);
  }
  if (strcmp(particle,"eta") == 0 && strcmp(sample,"data") == 0) {
    if (pt >= 6) h->Rebin(2);
    if (pt >= 8) h->Rebin(2);
  }
  RooHist * hpull;
  double *ypull;
  int numFitPar = 0;
  int ndf = 0;
  double chisq = 0;
  int nFullBinsPull = 0;
  int startbin;
  int endbin;
  int nBins;
  const char * sigtype;
  RooFitResult * res;
  RooGenericPdf * sig;
  RooGenericPdf * bkg;
  // Fitting with RooFit
  RooDataHist data("binnedData", "Binned dataset", RooArgSet(*(ws->var("mass"))), h);
  RooPlot * plot;

  int ntrials = 1;
  if (strcmp(particle,"pi0")==0) ntrials = 3;
  for (int trial = 0; trial < ntrials; trial++) {
    numFitPar = 0;
    ndf = 0;
    chisq = 0;
    nFullBinsPull = 0;
    if (trial > 0 ) { initialmean *= 1.03;  initialsigma *= 0.5; }
    plot = ws->var("mass")->frame(nbins);
    data.plotOn(plot, RooFit::DataError(RooAbsData::Poisson),RooFit::Name("data"));

    
    // Gaussian/double gaussian
    RooRealVar mu("mean", "mean of gaussian",initialmean,initialmean-0.1,initialmean+0.1);
    RooRealVar sigma1("sigma1", "width of gaussian", initialsigma, 0.005, 0.15); 
    RooRealVar sigma2("sigma2", "width of gaussian", initialsigma*1.1, 0.001, 0.2);
    RooGaussian * signal1 = new RooGaussian("signal1", "gaussian PDF", *(ws->var("mass")), mu, sigma1);
    RooGaussian * signal2 = new RooGaussian("signal2", "gaussian PDF", *(ws->var("mass")), mu, sigma2);
    RooRealVar gaussFrac("gaussFrac", "fraction of Gaussians", 0.5, 0, 1);
    RooAddPdf * doubleGauss = new RooAddPdf("doubleGauss", "Double Gaussian", RooArgList(*signal1,*signal2), gaussFrac);

    // Crystal Ball
    RooRealVar x0("x0", "x0", 0, initialmean-0.05, initialmean+0.05);
    RooRealVar sigma("sigma0", "sigma", initialsigma, 0.001, 0.2);
    RooRealVar alpha("alpha", "alpha", 10, -100, 100); // Negative for a right-side tail
    RooRealVar n("n", "n", 1, 0.1, 10);
    RooCrystalBall * cb = new RooCrystalBall("crystalBall", "Crystal Ball PDF", *(ws->var("mass")), x0, sigma, alpha, n);

    // polynomial/chebyshev background
    RooRealVar A("a", "1st order coefficient",Astart,Alow,Ahigh);
    RooRealVar B("b", "2nd order coefficient",Bstart,Blow,Bhigh);
    RooRealVar C("c", "3rd order coefficient",-100,-1e9,1e9);
    RooRealVar D("d", "4th order coefficient",0,-1e9,1e9);
    RooPolynomial *line = new RooPolynomial( "line",  "first order polynomial",  *(ws->var("mass")), RooArgList(A));
    RooPolynomial *poly2 = new RooPolynomial("poly2", "second order polynomial", *(ws->var("mass")), RooArgList(A,B));
    RooPolynomial *poly3 = new RooPolynomial("poly3", "third order polynomial",  *(ws->var("mass")), RooArgList(A,B,C));
    RooPolynomial *poly4 = new RooPolynomial("poly4", "fourth order polynomial", *(ws->var("mass")), RooArgList(A,B,C,D));
    RooChebychev *cheb2 = new RooChebychev("cheb", "second order chebychev", *(ws->var("mass")), RooArgList(A,B));
    RooChebychev *cheb3 = new RooChebychev("cheb", "third order chebychev",  *(ws->var("mass")), RooArgList(A,B,C));


    // Different choices for background and signal
    const char * bkgtype = map_bkg[textmap[particle]][textmap[sample]][textmap[trigger]][pt];
    sigtype = map_sig[textmap[particle]][textmap[sample]][textmap[trigger]][pt];

    if (strcmp("gauss",sigtype) == 0)       sig = (RooGenericPdf*) signal1;
    if (strcmp("doublegauss",sigtype) == 0) sig = (RooGenericPdf*) doubleGauss;
    if (strcmp("cb",sigtype) == 0)          sig = (RooGenericPdf*) cb;

    if (strcmp("line" ,bkgtype) == 0) bkg = (RooGenericPdf*) line;
    if (strcmp("poly2",bkgtype) == 0) bkg = (RooGenericPdf*) poly2;
    if (strcmp("poly3",bkgtype) == 0) bkg = (RooGenericPdf*) poly3;
    if (strcmp("poly4",bkgtype) == 0) bkg = (RooGenericPdf*) poly4;

    RooRealVar nsig("nsig", "signal count", 1000, 10, 1e8);
    RooRealVar nbkg("nbkg", "background count", 100, 0, 1e8);

    RooAddPdf * model = new RooAddPdf("model", "signal + background", RooArgList(*sig,*bkg),RooArgList(nsig, nbkg));
    ws->import(*model);

    RooLinkedList* fitcmd = new RooLinkedList();
    RooCmdArg opt1 = RooFit::Save(); fitcmd->Add(&opt1);
    //RooCmdArg opt2 = RooFit::PrefitDataFraction(0.1); fitcmd->Add(&opt2);
    RooCmdArg opt3 = RooFit::Minimizer("Minuit","minimize"); fitcmd->Add(&opt3);//"migradimproved");//);
    RooCmdArg opt4 = RooFit::NumCPU(16, 0); fitcmd->Add(&opt4);
    //RooCmdArg opt5 = RooFit::Range(h[i]->GetXaxis()->GetBinLowEdge(1), h[i]->GetXaxis()->GetBinUpEdge(h[i]->GetNbinsX())); fitcmd->Add(&opt5);
    RooCmdArg opt5 = RooFit::Range(plotlow,plothigh); fitcmd->Add(&opt5);
    RooCmdArg opt7 = RooFit::SumW2Error(kTRUE); fitcmd->Add(&opt7);
    RooCmdArg opt8 = RooFit::Extended(kTRUE); fitcmd->Add(&opt8);
    RooCmdArg opt9 = RooFit::Minos(kFALSE); fitcmd->Add(&opt9);
    RooCmdArg opt10= RooFit::Hesse(kTRUE); fitcmd->Add(&opt10);//RooCmdArg::none();
    //RooCmdArg opt11= RooFit::EvalErrorWall(kTRUE); fitcmd->Add(&opt11);//tRooCmdArg::none();
    //RooCmdArg opt12= RooFit::Strategy(2); fitcmd->Add(&opt12);//tRooCmdArg::none();
    //RooCmdArg opt13= RooFit::RecoverFromUndefinedRegions(3); fitcmd->Add(&opt13);//tRooCmdArg::none();
    //RooCmdArg opt14= RooFit::BatchMode(kTRUE); fitcmd->Add(&opt14);//tRooCmdArg::none();
    RooCmdArg opt15= RooFit::PrintEvalErrors(-1); fitcmd->Add(&opt15);//tRooCmdArg::none();

    res = ws->pdf("model")->fitTo(data, *fitcmd);
    ws->pdf("model")->plotOn(plot,Name("model"),LineColor(kRed),LineWidth(2),LineStyle(kSolid));
    ws->pdf("model")->plotOn(plot,Name("Sig"),Components(RooArgSet(*sig)),LineColor(kGreen+2),LineWidth(2),LineStyle(kSolid),FillStyle(3001),FillColor(kGreen+1),DrawOption("LF"));
    ws->pdf("model")->plotOn(plot,Name("Bkg"),Components(RooArgSet(*bkg)),LineColor(kBlue+2),LineStyle(kDashed),LineWidth(2));

    hpull = plot->pullHist("data","model");
    ypull = hpull->GetY();
    startbin = hpull->GetXaxis()->FindBin(plotlow);
    endbin = hpull->GetXaxis()->FindBin(plothigh);
    nBins = hpull->GetN(); 
    cout << "before-- chisq: " << chisq << " ndf: " << ndf << " chisq/ndf: " << chisq/ndf << endl;
    for(int i = 0; i < nBins; i++) {
      if (ypull[i] == 0 || std::isnan(ypull[i])) continue;
      chisq += TMath::Power(ypull[i],2);
      nFullBinsPull++;
    }
    numFitPar = res->floatParsFinal().getSize();
    ndf = nFullBinsPull - numFitPar;
    cout << "after-- chisq: " << chisq << " ndf: " << ndf << " chisq/ndf: " << chisq/ndf << endl;

    if (chisq/ndf < 8 || trial == 2) break;
    plot->remove("model");
    plot->remove("Sig");
    plot->remove("Bkg");
  }

  ws->import(*res);
  TPad * pplot = new TPad("pplot","",0,0.33,1,1);
  //TPad * pplot = new TPad("pplot","",0,0,1,1);
  pplot->SetBottomMargin(0.13);
  pplot->SetBottomMargin(0);
  pplot->SetRightMargin(0.02);
  pplot->SetLeftMargin(0.18);
  pplot->SetTopMargin(0.02);
  TPad * ppull = new TPad("ppull","",0,0,1,0.33);
  ppull->SetTopMargin(0);
  ppull->SetRightMargin(0.02);
  ppull->SetLeftMargin(0.18);
  ppull->SetBottomMargin(0.23);
  pplot->cd();
  gPad->SetTicks(1,1);

  double maxYValueFit = 0;

  for (int j = 0; j < plot->numItems(); j++) {
    if (plot->getObject(j)->InheritsFrom("RooCurve")) {
      RooCurve* curve = (RooCurve*)plot->getObject(j);
      if (strcmp(plot->getObject(j)->GetName(), "model") == 0) { 
        maxYValueFit = curve->getYAxisMax(); 
        break;
      }
    }
  }


  plot->SetFillStyle(4000);
  plot->GetXaxis()->SetLimits(plotlow, plothigh);
  plot->GetYaxis()->SetLimits(0,maxYValueFit*3);
  plot->GetYaxis()->SetRangeUser(0,maxYValueFit*3);
  plot->GetXaxis()->SetRangeUser(plotlow, plothigh);
  plot->GetYaxis()->SetTitleOffset(1.4);
  plot->GetYaxis()->CenterTitle();
  plot->GetYaxis()->SetTitleSize(0.048);
  plot->GetYaxis()->SetLabelSize(0.04);
  plot->GetXaxis()->SetLabelSize(0);
  plot->GetXaxis()->SetTitleSize(0);
  plot->SetTitle("");
  //plot->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");
  //plot->GetXaxis()->SetTitleOffset(1.20) ;
  //plot->GetXaxis()->SetLabelOffset(0.04) ;
  //plot->GetXaxis()->SetTitleSize(0.048) ;
  //plot->GetXaxis()->SetLabelSize(0.04) ;
  //plot->GetXaxis()->CenterTitle();
  plot->Draw();

  drawText("#bf{#it{sPHENIX}} Internal",0.55,0.81,1,22);
  const char * runtext;
  if (strcmp(sample,"data") == 0) runtext = Form("Run 47289-53879 %s",trigger);
  else if (strcmp(sample,"MC") == 0) runtext = Form("MC Pythia run 28 %s",trigger);
  drawText(runtext,0.55,0.74,1,22);
  drawText(Form("%i GeV < p_{T,#gamma#gamma} < %i GeV",pt,pt+1),0.6,0.68,1,16);
  drawText("|vz| < 30 cm",0.6,0.63,1,16);
  drawText("|#eta| < 1",0.6,0.58,1,16);
  drawText("p_{T,#gamma} > 1 GeV",0.6,0.53,1,16);
  drawText("#alpha < 0.6",0.6,0.48,1,16);

  ppull->cd();
  gPad->SetTicks(1,1);
  hpull->SetMarkerSize(0.7);
  RooPlot* pullFrame = ws->var("mass")->frame(Title(" ")) ;
  pullFrame->addPlotable(hpull,"P") ;
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.6) ;
  pullFrame->GetYaxis()->SetTitle("Pull") ;
  pullFrame->GetYaxis()->SetTitleSize(0.088) ;
  pullFrame->GetYaxis()->SetLabelSize(0.08) ;
  pullFrame->GetYaxis()->SetRangeUser(-16,16) ;
  pullFrame->GetXaxis()->SetLimits(plotlow, plothigh) ;
  pullFrame->GetXaxis()->SetRangeUser(plotlow, plothigh) ;
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");
  pullFrame->GetXaxis()->SetTitleOffset(1.20) ;
  pullFrame->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame->GetXaxis()->SetTitleSize(0.088) ;
  pullFrame->GetXaxis()->SetLabelSize(0.08) ;
  pullFrame->GetXaxis()->CenterTitle();

  pullFrame->GetYaxis()->SetTickSize(0.02);
  pullFrame->GetYaxis()->SetNdivisions(505);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw();

  if (startbin < 0 || endbin >= startbin + nBins) {
    std::cout << "startbin / endbin / nBins : " << startbin << ", " << endbin << ", " << nBins << std::endl;
    std::cerr << "Error: startbin or endbin is out of range." << std::endl;
    //return;
  }


  TLine *l1 = new TLine(plotlow,0,plothigh,0);
  l1->SetLineStyle(9);
  l1->Draw("same");



  pplot->cd();
  const char * muname;
  const char * stdname;
  if (strcmp("cb",sigtype) == 0) {
    muname = "x0";
    stdname = "sigma0";
  }
  else if (strcmp("gauss",sigtype) == 0 || strcmp("doublegauss",sigtype) == 0) {
    muname = "mean";
    stdname = "sigma1";
  }
  RooRealVar* nsigVar = (RooRealVar*) res->floatParsFinal().find("nsig");
  RooRealVar* muVar = (RooRealVar*) res->floatParsFinal().find(muname);
  RooRealVar* stdVar = (RooRealVar*) res->floatParsFinal().find(stdname);
  float nsigyield = nsigVar->getVal();
  float nsigyielderr = nsigVar->getError(); 
  float mean = muVar->getVal();
  float meanerr = muVar->getError();
  float std = stdVar->getVal();
  float stderr = stdVar->getError();

  drawText(Form("Fit #chi^2/NDF = %.2f",chisq/ndf),0.2,0.68,1,16);
  drawText(Form("#mu = %.4f#pm%.4f",mean,meanerr),0.2,0.63,1,16);
  drawText(Form("#sigma = %.4f#pm%.4f",std,stderr),0.2,0.58,1,16);
  drawText(Form("Yield: %.2f#pm%.2f",nsigyield,nsigyielderr),0.2,0.53,1,16);
  drawText(bkg->GetTitle(),0.2,0.48,1,16);
  c->cd();
  pplot->Draw();
  ppull->Draw();


  //return;
  c->SaveAs(Form("pdfs/mass_fits_%s_%s_pt%i_%s.pdf",sample,particle,pt,trigger));
  wf->cd();
  h->Write();
  wf->cd();
  ws->Write(); 
  workspacename = Form("workspaces/workspace_fits_%s_%s_pt%i_%s.root",sample,particle,pt,trigger);
  workspacesufficient = Form("sufficient_without_prob_1GeV_smear/workspace_fits_%s_%s_pt%i_%s.root",sample,particle,pt,trigger);
  pdffilename = Form("pdfs/mass_fits_%s_%s_pt%i_%s.pdf",sample,particle,pt,trigger);
  pdfsufficient = Form("sufficient_without_prob_1GeV_smear/mass_fits_%s_%s_pt%i_%s.pdf",sample,particle,pt,trigger);
  move();

}

