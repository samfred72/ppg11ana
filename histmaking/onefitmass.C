#include "commonUtility.h"
#include "Style_jaebeom.h"
#include <map>

using namespace RooFit;

const char * workspacename;
const char * workspacesufficient;
const char * pdffilename;
const char * pdfsufficient;

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


void onefitmass(const char * particle = "pi0", const char * type = "MB", const char *  subytpe = "", int pt = 1) {
  map<string,map<int,float>> map_plotlow;
  map<string,map<int,float>> map_plothigh;
  map<string,map<int,float>> map_initialmean;
  map<string,map<int,float>> map_initialsigma;
  map<string,map<int,float>> map_Astart;
  map<string,map<int,float>> map_Alow;
  map<string,map<int,float>> map_Ahigh;
  map_plotlow["pi0"][1] = 0.05; map_plotlow["eta"][1] = 0.35;
  map_plotlow["pi0"][2] = 0.05; map_plotlow["eta"][2] = 0.35;
  map_plotlow["pi0"][3] = 0.05; map_plotlow["eta"][3] = 0.35;
  map_plotlow["pi0"][4] = 0.05; map_plotlow["eta"][4] = 0.35;
  map_plotlow["pi0"][5] = 0.09; map_plotlow["eta"][5] = 0.35;
  map_plotlow["pi0"][6] = 0.09; map_plotlow["eta"][6] = 0.35;
  map_plothigh["pi0"][1] = 0.23; map_plothigh["eta"][1] = 0.75;
  map_plothigh["pi0"][2] = 0.23; map_plothigh["eta"][2] = 0.75;
  map_plothigh["pi0"][3] = 0.23; map_plothigh["eta"][3] = 0.75;
  map_plothigh["pi0"][4] = 0.23; map_plothigh["eta"][4] = 0.75;
  map_plothigh["pi0"][5] = 0.20; map_plothigh["eta"][5] = 0.75;
  map_plothigh["pi0"][6] = 0.20; map_plothigh["eta"][6] = 0.75;
  map_initialmean["pi0"][1] = 0.13; map_initialmean["eta"][1] = 0.60;
  map_initialmean["pi0"][2] = 0.13; map_initialmean["eta"][2] = 0.60;
  map_initialmean["pi0"][3] = 0.13; map_initialmean["eta"][3] = 0.60;
  map_initialmean["pi0"][4] = 0.13; map_initialmean["eta"][4] = 0.60;
  map_initialmean["pi0"][5] = 0.13; map_initialmean["eta"][5] = 0.60;
  map_initialmean["pi0"][6] = 0.13; map_initialmean["eta"][6] = 0.60;
  map_initialsigma["pi0"][1] = 0.01; map_initialsigma["eta"][1] = 0.03;
  map_initialsigma["pi0"][2] = 0.01; map_initialsigma["eta"][2] = 0.03;
  map_initialsigma["pi0"][3] = 0.01; map_initialsigma["eta"][3] = 0.03;
  map_initialsigma["pi0"][4] = 0.01; map_initialsigma["eta"][4] = 0.03;
  map_initialsigma["pi0"][5] = 0.01; map_initialsigma["eta"][5] = 0.03;
  map_initialsigma["pi0"][6] = 0.01; map_initialsigma["eta"][6] = 0.03;
  map_Astart["pi0"][1] = -10; map_Astart["eta"][1] = -10;
  map_Astart["pi0"][2] = -10; map_Astart["eta"][2] = -10;
  map_Astart["pi0"][3] = -10; map_Astart["eta"][3] = -10;
  map_Astart["pi0"][4] = -10; map_Astart["eta"][4] = -10;
  map_Astart["pi0"][5] = 0; map_Astart["eta"][5] = -10;
  map_Astart["pi0"][6] = 0; map_Astart["eta"][6] = 10;
  map_Alow["pi0"][1] = -1e7; map_Alow["eta"][1] = -1e1;
  map_Alow["pi0"][2] = -1e7; map_Alow["eta"][2] = -1e1;
  map_Alow["pi0"][3] = -1e7; map_Alow["eta"][3] = -1e1;
  map_Alow["pi0"][4] = -1e7; map_Alow["eta"][4] = -1e1;
  map_Alow["pi0"][5] = -1e7; map_Alow["eta"][5] = -1e1;
  map_Alow["pi0"][6] = -1e7; map_Alow["eta"][6] = -1e1;
  map_Ahigh["pi0"][1] = 1e7; map_Ahigh["eta"][1] = 1e1;
  map_Ahigh["pi0"][2] = 1e7; map_Ahigh["eta"][2] = 1e1;
  map_Ahigh["pi0"][3] = 1e7; map_Ahigh["eta"][3] = 1e1;
  map_Ahigh["pi0"][4] = 1e7; map_Ahigh["eta"][4] = 1e1;
  map_Ahigh["pi0"][5] = 1e7; map_Ahigh["eta"][5] = 1e1;
  map_Ahigh["pi0"][6] = 1e7; map_Ahigh["eta"][6] = 1e1;

  float plotlow =      map_plotlow[particle][pt];
  float plothigh =     map_plothigh[particle][pt];
  float initialmean =  map_initialmean[particle][pt];
  float initialsigma = map_initialsigma[particle][pt];
  float Alow =         map_Alow[particle][pt];
  float Ahigh =        map_Ahigh[particle][pt];
  float Astart =       map_Astart[particle][pt];

  TFile * wf = new TFile(Form("workspaces/workspace_fits_%s_%s_pt%i.root",type,particle,pt),"RECREATE"); 
  TFile * inf = TFile::Open(Form("hists/mass_%s.root",type),"READ");
  TH1D * h = shorthist((TH1D*)inf->Get(Form("hmass_%s_pt%i",particle,pt)),plotlow,plothigh);

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

  if (strcmp(particle,"eta") == 0 && strcmp(type,"MB") == 0) {
    h->Rebin(2);
    if (pt > 4) h->Rebin(2);
  }
  if (strcmp(particle,"pi0") == 0 && strcmp(type,"MB") == 0) {
    if (pt > 5) h->Rebin(2);
  }
  // Fitting with RooFit
  RooDataHist data("binnedData", "Binned dataset", RooArgSet(*(ws->var("mass"))), h);
  RooPlot * plot = ws->var("mass")->frame(nbins);
  data.plotOn(plot, RooFit::DataError(RooAbsData::Poisson),RooFit::Name("data"));

  // Gaussian/double gaussian
  RooRealVar mu("mean", "mean of gaussian",initialmean,initialmean-0.05,initialmean+0.05);
  RooRealVar sigma1("sigma1", "width of gaussian", initialsigma, 0.005, 0.07); 
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
  RooRealVar A("a", "1th order coefficient",Astart,Alow,Ahigh);
  RooRealVar B("b", "2th order coefficient",-10,-1e8,1e8);
  RooRealVar C("c", "3th order coefficient",-100,-1e9,1e9);
  RooRealVar D("d", "4th order coefficient",0,-1e9,1e9);
  RooPolynomial *line = new RooPolynomial("line", "first order polynomial", *(ws->var("mass")), RooArgList(A));
  RooPolynomial *poly2 = new RooPolynomial("poly2", "second order polynomial", *(ws->var("mass")), RooArgList(A,B));
  RooPolynomial *poly3 = new RooPolynomial("poly3", "third order polynomial", *(ws->var("mass")), RooArgList(A,B,C));
  RooPolynomial *poly4 = new RooPolynomial("poly4", "fourth order polynomial", *(ws->var("mass")), RooArgList(A,B,C,D));
  RooChebychev *cheb2 = new RooChebychev("cheb", "second order chebychev", *(ws->var("mass")), RooArgList(A,B));
  RooChebychev *cheb3 = new RooChebychev("cheb", "third order chebychev", *(ws->var("mass")), RooArgList(A,B,C));

  RooGenericPdf * bkg;
  // Different choices for background and signal
  if (strcmp(particle,"pi0") == 0) {
    if (strcmp(type,"MB") == 0) bkg = (RooGenericPdf*) poly2;
    if (strcmp(type,"Jet10") == 0) bkg = (RooGenericPdf*) poly3;
    if (strcmp(type,"Jet20") == 0) bkg = (RooGenericPdf*) poly3;
    if (strcmp(type,"Jet30") == 0) bkg = (RooGenericPdf*) poly3;
  }
  else {
    if (strcmp(type,"MB") == 0) bkg = (RooGenericPdf*) line;
    if (strcmp(type,"Jet10") == 0) bkg = (RooGenericPdf*) poly2;
    if (strcmp(type,"Jet20") == 0) bkg = (RooGenericPdf*) poly2;
    if (strcmp(type,"Jet30") == 0) bkg = (RooGenericPdf*) poly2;
  }

  RooGenericPdf * sig;
  const char * sigtype = "gauss";
  if (strcmp("gauss",sigtype) == 0)       sig = (RooGenericPdf*) signal1;
  if (strcmp("doublegauss",sigtype) == 0) sig = (RooGenericPdf*) doubleGauss;
  if (strcmp("cb",sigtype) == 0)          sig = (RooGenericPdf*) cb;

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
  RooCmdArg opt14= RooFit::BatchMode(kTRUE); fitcmd->Add(&opt14);//tRooCmdArg::none();
  RooCmdArg opt15= RooFit::PrintEvalErrors(-1); fitcmd->Add(&opt15);//tRooCmdArg::none();

  RooFitResult * res = ws->pdf("model")->fitTo(data, *fitcmd);
  ws->import(*res);



  ws->pdf("model")->plotOn(plot,Name("model"),LineColor(kRed),LineWidth(2),LineStyle(kSolid));
  ws->pdf("model")->plotOn(plot,Name("Sig"),Components(RooArgSet(*sig)),LineColor(kGreen+2),LineWidth(2),LineStyle(kSolid),FillStyle(3001),FillColor(kGreen+1),DrawOption("LF"));
  ws->pdf("model")->plotOn(plot,Name("Bkg"),Components(RooArgSet(*bkg)),LineColor(kBlue+2),LineStyle(kDashed),LineWidth(2));


  TPad * pplot = new TPad("pplot","",0,0.33,1,1);
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
  plot->Draw();

  drawText("#bf{#it{sPHENIX}} Internal",0.55,0.81,1,22);
  const char * runtext = Form("MC Pythia run 28 %s",type);
  drawText(runtext,0.55,0.74,1,22);
  drawText(Form("%i GeV < p_{T,#gamma#gamma} < %i GeV",pt,pt+1),0.6,0.68,1,16);
  drawText("|vz| < 30 cm",0.6,0.63,1,16);
  drawText("|#eta| < 1",0.6,0.58,1,16);
  drawText("p_{T,#gamma} > 0.5 GeV",0.6,0.53,1,16);
  drawText("prob_{#gamma} > 0.05",0.6,0.48,1,16);
  drawText("#alpha < 0.6",0.6,0.43,1,16);

  ppull->cd();
  gPad->SetTicks(1,1);
  RooHist * hpull = plot->pullHist("data","model");
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

  double chisq = 0;
  int nFullBinsPull = 0;
  int startbin = hpull->GetXaxis()->FindBin(plotlow);
  int endbin = hpull->GetXaxis()->FindBin(plothigh);
  int nBins = hpull->GetN(); 
  double *ypull = hpull->GetY();
  if (startbin < 0 || endbin >= startbin + nBins) {
    std::cout << "startbin / endbin / nBins : " << startbin << ", " << endbin << ", " << nBins << std::endl;
    std::cerr << "Error: startbin or endbin is out of range." << std::endl;
    //return;
  }

  for(int i = 0; i < nBins; i++) {
    if (ypull[i] == 0 || std::isnan(ypull[i])) continue;
    chisq += TMath::Power(ypull[i],2);
    nFullBinsPull++;
  }
  int numFitPar = res->floatParsFinal().getSize();
  int ndf = nFullBinsPull - numFitPar;

  TLine *l1 = new TLine(plotlow,0,plothigh,0);
  l1->SetLineStyle(9);
  l1->Draw("same");

  cout << "chisq: " << chisq << " ndf: " << ndf << " chisq/ndf: " << chisq/ndf << endl;


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
  c->SaveAs(Form("pdfs/mass_fits_%s_%s_pt%i.pdf",type,particle,pt));
  wf->cd();
  h->Write();
  wf->cd();
  ws->Write(); 
  workspacename = Form("workspaces/workspace_fits_%s_%s_pt%i.root",type,particle,pt);
  workspacesufficient = Form("sufficient/workspace_fits_%s_%s_pt%i.root",type,particle,pt);
  pdffilename = Form("pdfs/mass_fits_%s_%s_pt%i.pdf",type,particle,pt);
  pdfsufficient = Form("sufficient/mass_fits_%s_%s_pt%i.pdf",type,particle,pt);

}

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
