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


void alphafitmass(const char * particle = "pi0", const char * sample = "MC", const char * trigger = "MB") {
  set_maps();

  TCanvas * c = new TCanvas("c","",600,800);
  gStyle->SetOptStat(0);
  c->SaveAs("alphafits.pdf[");
  map<string,int> textmap = {
    {"pi0",0},{"eta",1},
    {"data",0},{"MC",1},
    {"MB",0},{"photon3",1},{"photon4",2},{"photon5",3},
    {"Jet10",1},{"Jet20",2},{"Jet30",3}
  };

  float ptbins[4] = {2,2.5,3,20};
  float abins[4] = {0,0.08,0.2,1};
  for (int pt = 0; pt < 3; pt++) { 
    for (int alph = 0; alph < 3; alph++) {
      for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
          for (int k = 0; k < 5; k++) {
            for (int l = 0; l < 5; l++) {
              float plotlow =      0.07;//map_plotlow     [textmap[particle]][textmap[sample]][textmap[trigger]][2];
              float plothigh =     map_plothigh    [textmap[particle]][textmap[sample]][textmap[trigger]][2];
              float initialmean =  .13;//map_initialmean =1.3;//[textmap[particle]][textmap[sample]][textmap[trigger]][pt];
              float initialsigma = 0.05;//map_initialsigma=0.05;//[textmap[particle]][textmap[sample]][textmap[trigger]][pt];
              float Alow =         map_Alow        [textmap[particle]][textmap[sample]][textmap[trigger]][2];
              float Ahigh =        map_Ahigh       [textmap[particle]][textmap[sample]][textmap[trigger]][2];
              float Astart =       map_Astart      [textmap[particle]][textmap[sample]][textmap[trigger]][2];
              float Blow =         map_Blow        [textmap[particle]][textmap[sample]][textmap[trigger]][2];
              float Bhigh =        map_Bhigh       [textmap[particle]][textmap[sample]][textmap[trigger]][2];
              float Bstart =       map_Bstart      [textmap[particle]][textmap[sample]][textmap[trigger]][2];

              TFile * inf;
              if (strcmp(sample,"data") == 0) inf = TFile::Open("data_hists/mass_data.root","READ");
              else if (strcmp(sample,"MC") == 0) inf = TFile::Open(Form("hists/mass_%s_%s.root",sample,trigger),"READ");
              TH1D * h;
              if (strcmp(sample,"data") == 0) h = shorthist((TH1D*)inf->Get(Form("halpha%i_%i",pt,alph)),plotlow,plothigh);
              else if (strcmp(sample,"MC") == 0) h = shorthist((TH1D*)inf->Get(Form("halpha%i_%i",pt,alph)),plotlow,plothigh);

              TFile * wf = new TFile(Form("truthfits_smeartest/workspace_pt%i_a%i_%i_%i_%i_%i.root",pt,alph,i,j,k,l),"RECREATE"); 
              RooWorkspace * ws;
              for (int trial = 0; trial < 3; trial++) {
                cout << "attempting trial: " << trial << endl;
                if (trial > 0 ) { initialmean *= 1.04;  initialsigma *= 0.5; }
                int binlow = h->FindBin(plotlow);
                int binhigh = h->FindBin(plothigh);
                int nbins = binhigh - binlow;
                cout << "Low bin: " << binlow << " High bin: " << binhigh << " nbins: " << nbins << endl;
                RooMsgService::instance().setSilentMode(true);
                ws = new RooWorkspace("workspace");
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
                // Fitting with RooFit
                RooDataHist data("binnedData", "Binned dataset", RooArgSet(*(ws->var("mass"))), h);
                RooPlot * plot = ws->var("mass")->frame(nbins);
                data.plotOn(plot, RooFit::DataError(RooAbsData::Poisson),RooFit::Name("data"));

                // Gaussian/double gaussian
                RooRealVar mu("mean", "mean of gaussian",initialmean,initialmean-0.05,initialmean+0.05);
                RooRealVar sigma1("sigma1", "width of gaussian", initialsigma, 0.005, 0.1); 
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
                RooGenericPdf * sig;
                RooGenericPdf * bkg;
                const char * bkgtype = map_bkg[textmap[particle]][textmap[sample]][textmap[trigger]][pt];
                const char * sigtype = map_sig[textmap[particle]][textmap[sample]][textmap[trigger]][pt];

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

                RooFitResult * res = ws->pdf("model")->fitTo(data, *fitcmd);
                ws->pdf("model")->plotOn(plot,Name("model"),LineColor(kRed),LineWidth(2),LineStyle(kSolid));
                ws->pdf("model")->plotOn(plot,Name("Sig"),Components(RooArgSet(*sig)),LineColor(kGreen+2),LineWidth(2),LineStyle(kSolid),FillStyle(3001),FillColor(kGreen+1),DrawOption("LF"));
                ws->pdf("model")->plotOn(plot,Name("Bkg"),Components(RooArgSet(*bkg)),LineColor(kBlue+2),LineStyle(kDashed),LineWidth(2));

                TPad * pplot = new TPad("pplot","",0,0,1,1);
                pplot->SetBottomMargin(0.13);
                pplot->SetRightMargin(0.02);
                pplot->SetLeftMargin(0.18);
                pplot->SetTopMargin(0.02);
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
                plot->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");
                plot->GetXaxis()->SetTitleOffset(1.20) ;
                plot->GetXaxis()->SetLabelOffset(0.04) ;
                plot->GetXaxis()->SetTitleSize(0.048) ;
                plot->GetXaxis()->SetLabelSize(0.04) ;
                plot->GetXaxis()->CenterTitle();
                plot->Draw();
                drawText("#bf{#it{sPHENIX}} Internal",0.55,0.81,1,22);
                const char * runtext;
                if (strcmp(sample,"data") == 0) runtext = Form("Run 47289-53879 %s",trigger);
                else if (strcmp(sample,"MC") == 0) runtext = Form("MC Pythia run 28 %s",trigger);
                drawText(runtext,0.55,0.74,1,22);
                drawText(Form("%f GeV < p_{T,#gamma#gamma} < %f GeV",ptbins[pt],ptbins[pt+1]),0.6,0.68,1,16);
                drawText("|vz| < 30 cm",0.6,0.63,1,16);
                drawText("|#eta| < 0.3",0.6,0.58,1,16);
                drawText("p_{T,#gamma} > 1 GeV",0.6,0.53,1,16);
                drawText("#alpha < 0.6",0.6,0.48,1,16);

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
                RooHist * hpull = plot->pullHist("data","model");
                double *ypull = hpull->GetY();
                int nFullBinsPull = 0;
                float chisq = 0;
                for(int i = 0; i < nBins; i++) {
                  if (ypull[i] == 0 || std::isnan(ypull[i])) continue;
                  chisq += TMath::Power(ypull[i],2);
                  nFullBinsPull++;
                }
                int numFitPar = res->floatParsFinal().getSize();
                int ndf = nFullBinsPull - numFitPar;
                drawText(Form("chisq/ndf: %.2f",chisq/ndf),0.2,0.73,1,16);
                drawText(Form("%f < #alpha < %f",abins[alph],abins[alph+1]),0.2,0.68,1,16);
                drawText(Form("#mu = %.4f#pm%.4f",mean,meanerr),0.2,0.63,1,16);
                drawText(Form("#sigma = %.4f#pm%.4f",std,stderr),0.2,0.58,1,16);
                drawText(Form("Yield: %.2f#pm%.2f",nsigyield,nsigyielderr),0.2,0.53,1,16);
                drawText(bkg->GetTitle(),0.2,0.48,1,16);
                c->cd();
                pplot->Draw();

                if (chisq/ndf < 10 && chisq/ndf > 0.05) break;
              }
              c->SaveAs("alphafits.pdf");
              ws->Write();
              delete ws;
              wf->Close();
              delete wf;
            }
          }
        }
      }
    }
  }
  c->SaveAs("alphafits.pdf]");
}

