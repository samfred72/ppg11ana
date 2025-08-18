#include "commonUtility.h"
#include "Style_jaebeom.h"

using namespace RooFit;

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


void fitmass(const char * particle = "pi0", const char * type = "MB",int starti = 0, int numplots = -1) {
  float plotlow =      (strcmp(particle,"pi0") == 0) ? 0.06 : 0.37;
  float plothigh =     (strcmp(particle,"pi0") == 0) ? 0.23 : 0.73;
  float initialmean =  (strcmp(particle,"pi0") == 0) ? 0.13 : 0.55;
  float initialsigma = (strcmp(particle,"pi0") == 0) ? 0.013 : 0.03;
  float Alow =         (strcmp(particle,"pi0") == 0) ? -1e9 : -1e2;
  float Ahigh =        (strcmp(particle,"pi0") == 0) ? 1e9 : 1e2;
  
  TFile * wf = new TFile(Form("workspace_fits_%s_%s.root",type,particle),"RECREATE"); 
  TFile * inf = TFile::Open(Form("mass_%s.root",type),"READ");
  const static int nhists = 25;
  TH1D * h[nhists];
  for (int i = 0; i < nhists; i++) {
    h[i] = shorthist((TH1D*)inf->Get(Form("hmass_%s_pt%i",particle,i)),plotlow,plothigh);
  }

  TH1D * p[nhists];
  TCanvas * c = new TCanvas("c","",600,800);
  c->SaveAs(Form("mass_fits_%s_%s.pdf[",type,particle));
  TPad * pplot[nhists];
  TPad * ppull[nhists];
  TLegend * leg[nhists];
  gStyle->SetOptStat(0);
  RooGenericPdf * sig[nhists];
  RooGenericPdf * bkg[nhists];
  RooAddPdf * model[nhists];
  RooPlot * plot[nhists];

  int binlow = h[0]->FindBin(plotlow);
  int binhigh = h[0]->FindBin(plothigh);
  int nbins = binhigh - binlow;
  cout << "Low bin: " << binlow << " High bin: " << binhigh << " nbins: " << nbins << endl;
  //RooMsgService::instance().setSilentMode(true);
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->factory("mass[0.0, 1.0]");
  ws->var("mass")->setRange(plotlow,plothigh);
  ws->var("mass")->Print();
  ws->var("mass")->setBins(nbins);
 
  if (numplots == -1) numplots = nhists;
  if (starti+numplots > nhists) numplots = nhists - starti;
  for (int i = starti; i < starti + numplots; i++) { 
    float integral = h[i]->Integral(1,h[i]->GetNbinsX());
    if (integral < 100) continue;
    
    if (strcmp(particle,"eta") == 0 && strcmp(type,"MB")) h[i]->Rebin(2);
    // Fitting with RooFit
    RooDataHist data("binnedData", "Binned dataset", RooArgSet(*(ws->var("mass"))), h[i]);
    plot[i] = ws->var("mass")->frame(nbins);
    data.plotOn(plot[i], RooFit::DataError(RooAbsData::Poisson),RooFit::Name("data"));
    
    // Gaussian/double gaussian
    RooRealVar mu(Form("mean%i",i), "mean of gaussian",initialmean,initialmean-0.05,initialmean+0.05);
    RooRealVar sigma1(Form("sigma1_%i",i), "width of gaussian", initialsigma, 0.001, 0.2); 
    RooRealVar sigma2(Form("sigma2_%i",i), "width of gaussian", initialsigma*1.1, 0.001, 0.2);
    RooGaussian * signal1 = new RooGaussian(Form("signal1_%i",i), "gaussian PDF", *(ws->var("mass")), mu, sigma1);
    RooGaussian * signal2 = new RooGaussian(Form("signal2_%i",i), "gaussian PDF", *(ws->var("mass")), mu, sigma2);
    RooRealVar gaussFrac(Form("gaussFrac%i",i), "fraction of Gaussians", 0.5, 0, 1);
    RooAddPdf * doubleGauss = new RooAddPdf(Form("doubleGauss%i",i), "Double Gaussian", RooArgList(*signal1,*signal2), gaussFrac);

    // Crystal Ball
    RooRealVar x0(Form("x0_%i",i), "x0", 0, initialmean-0.05, initialmean+0.05);
    RooRealVar sigma(Form("sigma0_%i",i), "sigma", initialsigma, 0.001, 0.2);
    RooRealVar alpha(Form("alpha%i",i), "alpha", 10, -100, 100); // Negative for a right-side tail
    RooRealVar n("n", "n", 1, 0.1, 10);
    RooCrystalBall * cb = new RooCrystalBall(Form("crystalBall%i",i), "Crystal Ball PDF", *(ws->var("mass")), x0, sigma, alpha, n);
    
    // polynomial/chebyshev background
    RooRealVar A(Form("a%i",i), "1th order coefficient",20000,Alow,Ahigh);
    RooRealVar B(Form("b%i",i), "2th order coefficient",-40000,-1e9,1e9);
    RooRealVar C(Form("c%i",i), "3th order coefficient",-100,-1e9,1e9);
    RooRealVar D(Form("d%i",i), "4th order coefficient",0,-1e9,1e9);
    RooPolynomial *poly2 = new RooPolynomial(Form("poly2_%i",i), "second order polynomial", *(ws->var("mass")), RooArgList(A,B));
    RooPolynomial *poly3 = new RooPolynomial(Form("poly3_%i",i), "third order polynomial", *(ws->var("mass")), RooArgList(A,B,C));
    RooPolynomial *poly4 = new RooPolynomial(Form("poly4_%i",i), "fourth order polynomial", *(ws->var("mass")), RooArgList(A,B,C,D));
    RooChebychev *cheb2 = new RooChebychev(Form("cheb_2%i",i), "second order chebychev", *(ws->var("mass")), RooArgList(A,B));
    RooChebychev *cheb3 = new RooChebychev(Form("cheb_3%i",i), "third order chebychev", *(ws->var("mass")), RooArgList(A,B,C));
    RooPolynomial *line = new RooPolynomial(Form("line%i",i), "first order polynomial", *(ws->var("mass")), RooArgList(A));


    // Different choices for background and signal
    if (strcmp(particle,"pi0") == 0) {
      if (strcmp(type,"MB") == 0) bkg[i] = (RooGenericPdf*) poly2;
      if (strcmp(type,"Jet10") == 0) bkg[i] = (RooGenericPdf*) poly3;
      if (strcmp(type,"Jet20") == 0) bkg[i] = (RooGenericPdf*) poly3;
      if (strcmp(type,"Jet30") == 0) bkg[i] = (RooGenericPdf*) poly3;
    }
    else {
      if (strcmp(type,"MB") == 0) bkg[i] = (RooGenericPdf*) line;
      if (strcmp(type,"Jet10") == 0) bkg[i] = (RooGenericPdf*) poly2;
      if (strcmp(type,"Jet20") == 0) bkg[i] = (RooGenericPdf*) poly2;
      if (strcmp(type,"Jet30") == 0) bkg[i] = (RooGenericPdf*) poly2;
    }
    const char * sigtype = "doublegauss";
    if (strcmp("gauss",sigtype) == 0)       sig[i] = (RooGenericPdf*) signal1;
    if (strcmp("doublegauss",sigtype) == 0) sig[i] = (RooGenericPdf*) doubleGauss;
    if (strcmp("cb",sigtype) == 0)          sig[i] = (RooGenericPdf*) cb;

    RooRealVar nsig(Form("nsig%i",i), "signal count", 1000, 10, 1e8);
    RooRealVar nbkg(Form("nbkg%i",i), "background count", 100, 0, 1e8);

    model[i] = new RooAddPdf(Form("model%i",i), "signal + background", RooArgList(*sig[i],*bkg[i]),RooArgList(nsig, nbkg));
    ws->import(*model[i]);

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
    
    cout << "on index " << i << endl;
    RooFitResult * res = ws->pdf(Form("model%i",i))->fitTo(data, *fitcmd);
    cout << "success fitting!" << endl;
    

    cout <<" HELLO???????????????????????????????????????????????" << endl;
    ws->pdf(Form("model%i",i))->plotOn(plot[i],Name(Form("model%i",i)),LineColor(kRed),LineWidth(2),LineStyle(kSolid));
    ws->pdf(Form("model%i",i))->plotOn(plot[i],Name(Form("Sig%i",i)),Components(RooArgSet(*sig[i])),LineColor(kGreen+2),LineWidth(2),LineStyle(kSolid),FillStyle(3001),FillColor(kGreen+1),DrawOption("LF"));
    ws->pdf(Form("model%i",i))->plotOn(plot[i],Name(Form("Bkg%i",i)),Components(RooArgSet(*bkg[i])),LineColor(kBlue+2),LineStyle(kDashed),LineWidth(2));
    

    pplot[i] = new TPad(Form("pplot%i",i),"",0,0.33,1,1);
    pplot[i]->SetBottomMargin(0);
    pplot[i]->SetRightMargin(0.02);
    pplot[i]->SetLeftMargin(0.18);
    pplot[i]->SetTopMargin(0.02);
    ppull[i] = new TPad(Form("ppull%i",i),"",0,0,1,0.33);
    ppull[i]->SetTopMargin(0);
    ppull[i]->SetRightMargin(0.02);
    ppull[i]->SetLeftMargin(0.18);
    ppull[i]->SetBottomMargin(0.23);
    pplot[i]->cd();
    gPad->SetTicks(1,1);
    
    double maxYValueFit = 0;

    for (int j = 0; j < plot[i]->numItems(); j++) {
      if (plot[i]->getObject(j)->InheritsFrom("RooCurve")) {
        RooCurve* curve = (RooCurve*)plot[i]->getObject(j);
        if (strcmp(plot[i]->getObject(j)->GetName(), Form("model%i",i)) == 0) { 
          maxYValueFit = curve->getYAxisMax(); 
          break;
        }
      }
    }


    plot[i]->SetFillStyle(4000);
    plot[i]->GetXaxis()->SetLimits(plotlow, plothigh);
    plot[i]->GetYaxis()->SetLimits(0,maxYValueFit*3);
    plot[i]->GetYaxis()->SetRangeUser(0,maxYValueFit*3);
    plot[i]->GetXaxis()->SetRangeUser(plotlow, plothigh);
    plot[i]->GetYaxis()->SetTitleOffset(1.4);
    plot[i]->GetYaxis()->CenterTitle();
    plot[i]->GetYaxis()->SetTitleSize(0.048);
    plot[i]->GetYaxis()->SetLabelSize(0.04);
    plot[i]->GetXaxis()->SetLabelSize(0);
    plot[i]->GetXaxis()->SetTitleSize(0);
    plot[i]->SetTitle("");
    plot[i]->Draw();
    
    drawText("#bf{#it{sPHENIX}} Internal",0.55,0.81,1,22);
    const char * runtext = (strcmp("MC",type) == 0) ? Form("MC Pythia run 21 MB") : Form("MC Pythia run 28 %s",type);
    drawText(runtext,0.55,0.74,1,22);
    drawText(Form("%i GeV < p_{T,#gamma#gamma} < %i GeV",i,i+1),0.6,0.68,1,16);
    drawText("|vz| < 30 cm",0.6,0.63,1,16);
    drawText("|#eta| < 1.1",0.6,0.58,1,16);
    drawText("p_{T,#gamma} > 0.5 GeV",0.6,0.53,1,16);
    drawText("prob_{#gamma} > 0.05",0.6,0.48,1,16);
    drawText("#alpha < 0.6",0.6,0.43,1,16);
    
    ppull[i]->cd();
    gPad->SetTicks(1,1);
    RooHist * hpull = plot[i]->pullHist("data",Form("model%i",i));
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


    pplot[i]->cd();
    const char * muname;
    const char * stdname;
    if (strcmp("cb",sigtype) == 0) {
      muname = Form("x0_%i",i);
      stdname = Form("sigma0_%i",i);
    }
    else if (strcmp("gauss",sigtype) == 0 || strcmp("doublegauss",sigtype) == 0) {
      muname = Form("mean%i",i);
      stdname = Form("sigma1_%i",i);
    }
    RooRealVar* nsigVar = (RooRealVar*) res->floatParsFinal().find(Form("nsig%i",i));
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
    c->cd();
    pplot[i]->Draw();
    ppull[i]->Draw();


    //return;
    c->SaveAs(Form("mass_fits_%s_%s.pdf",type,particle));
    wf->cd();
    h[i]->Write();
  }
  c->SaveAs(Form("mass_fits_%s_%s.pdf]",type,particle));
  wf->cd();
  ws->Write(); 
}
    
