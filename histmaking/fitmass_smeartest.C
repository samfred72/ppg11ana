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


void fitmass_smeartest(const char * particle = "pi0", const char * sample = "MC", const char * trigger = "MB") {
  set_maps();

  TFile * inf = TFile::Open(Form("mass_smeartest.root"),"READ");
  TCanvas * c = new TCanvas("c","",600,800);
  gStyle->SetOptStat(0);
  c->SaveAs("smearfits.pdf[");
  map<string,int> textmap = {
    {"pi0",0},{"eta",1},
    {"data",0},{"MC",1},
    {"MB",0},{"photon3",1},{"photon4",2},{"photon5",3},
    {"Jet10",1},{"Jet20",2},{"Jet30",3}
  };

  cout << "getting data" << endl;
  float datapeaks[3] = {0,0,0};
  float datawidths[3] = {0,0,0};
  /*
  for (int i = 0; i < 3; i++) {
    TFile * f;
    RooWorkspace *w;
    RooFitResult *fitRes;
    RooRealVar *mean;
    RooRealVar *sigma;
    RooFormulaVar *ratio;
    const char * filename = Form("/sphenix/user/samfred/run25/ppg11/histmaking/sufficient_without_prob_1GeV/workspace_fits_data_%s_pt%i_mbd.root",particle,i);
    if (!gSystem->AccessPathName(filename)) {
      f = TFile::Open(filename,"READ");
      w = (RooWorkspace*) f->Get("workspace");
      mean = w->var("mean"); 
      sigma = w->var("sigma1"); 
      datapeaks[i] = mean->getVal();
      datawidths[i] = sigma->getVal();
    }
    delete f;
  }

  */
  cout << "doing fits" << endl;
  float fitpeaks[3][3][5][5][5][5];
  float fitwidths[3][3][5][5][5][5];
  float fitchisq[3][3][5][5][5][5];
  for (int pt = 0; pt < 3; pt++) { 
    float plotlow =      map_plotlow     [textmap[particle]][textmap[sample]][textmap[trigger]][2];
    float plothigh =     map_plothigh    [textmap[particle]][textmap[sample]][textmap[trigger]][2];
    float initialmean =  0.135;//map_initialmean [textmap[particle]][textmap[sample]][textmap[trigger]][pt];
    float initialsigma = 0.05;//map_initialsigma[textmap[particle]][textmap[sample]][textmap[trigger]][pt];
    float Alow =         map_Alow        [textmap[particle]][textmap[sample]][textmap[trigger]][2];
    float Ahigh =        map_Ahigh       [textmap[particle]][textmap[sample]][textmap[trigger]][2];
    float Astart =       map_Astart      [textmap[particle]][textmap[sample]][textmap[trigger]][2];
    float Blow =         map_Blow        [textmap[particle]][textmap[sample]][textmap[trigger]][2];
    float Bhigh =        map_Bhigh       [textmap[particle]][textmap[sample]][textmap[trigger]][2];
    float Bstart =       map_Bstart      [textmap[particle]][textmap[sample]][textmap[trigger]][2];
    for (int a = 0; a < 3; a++) {
      for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
          for (int k = 0; k < 5; k++) {
            for (int l = 0; l < 5; l++) {
              TFile * wf = TFile::Open(Form("truthfits_smeartest/workspace_pt%i_a%i_%i_%i_%i_%i.root",pt,a,i,j,k,l),"RECREATE");
              initialmean =  map_initialmean [textmap[particle]][textmap[sample]][textmap[trigger]][pt];
              initialsigma = map_initialsigma[textmap[particle]][textmap[sample]][textmap[trigger]][pt];
              for (int trial = 0; trial < 3; trial++) {
                if (trial > 0 ) { initialmean *= 1.03;  initialsigma *= 0.5; }
                cout << "Fitting " << pt << " " << a << " " << i << " " << j << " " << k << " " << l << endl;
                TH1D * h = shorthist((TH1D*)inf->Get(Form("hmass_%s_pt%i_a%i_%i_%i_%i_%i",particle,pt,a,i,j,k,l)),plotlow,plothigh);

                int binlow = h->FindBin(plotlow);
                int binhigh = h->FindBin(plothigh);
                int nbins = binhigh - binlow;
                cout << "Low bin: " << binlow << " High bin: " << binhigh << " nbins: " << nbins << endl;
                RooMsgService::instance().setSilentMode(true);
                RooWorkspace * ws = new RooWorkspace("workspace");
                ws->factory("mass[0.0, 1.0]");
                ws->var("mass")->setRange(plotlow,plothigh);
                ws->var("mass")->Print();
                ws->var("mass")->setBins(nbins);

                float integral = h->Integral(1,h->GetNbinsX());
                if (integral < 10) {
                  cout << "Histogram too empty! Exiting." << endl;
                  trial = 3;
                  continue;
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
                RooGaussian * signal1 = new RooGaussian("signal1", "gaussian PDF", *(ws->var("mass")), mu, sigma1);

                // polynomial/chebyshev background
                RooRealVar A("a", "1st order coefficient",Astart,Alow,Ahigh);
                RooRealVar B("b", "2nd order coefficient",Bstart,Blow,Bhigh);
                RooRealVar C("c", "3rd order coefficient",-100,-1e9,1e9);
                RooRealVar D("d", "4th order coefficient",0,-1e9,1e9);
                RooPolynomial *line = new RooPolynomial( "line",  "first order polynomial",  *(ws->var("mass")), RooArgList(A));
                RooPolynomial *poly2 = new RooPolynomial("poly2", "second order polynomial", *(ws->var("mass")), RooArgList(A,B));
                RooPolynomial *poly3 = new RooPolynomial("poly3", "third order polynomial",  *(ws->var("mass")), RooArgList(A,B,C));

                // Different choices for background and signal
                RooGenericPdf * sig;
                RooGenericPdf * bkg;
                const char * bkgtype = map_bkg[textmap[particle]][textmap[sample]][textmap[trigger]][pt];

                sig = (RooGenericPdf*) signal1;

                if (strcmp("line" ,bkgtype) == 0) bkg = (RooGenericPdf*) line;
                if (strcmp("poly2",bkgtype) == 0) bkg = (RooGenericPdf*) poly2;
                if (strcmp("poly3",bkgtype) == 0) bkg = (RooGenericPdf*) poly3;

                RooRealVar nsig("nsig", "signal count", 1000, 10, 1e8);
                RooRealVar nbkg("nbkg", "background count", 100, 0, 1e8);

                RooAddPdf * model = new RooAddPdf("model", "signal + background", RooArgList(*sig,*bkg),RooArgList(nsig, nbkg));
                ws->import(*model);

                RooLinkedList* fitcmd = new RooLinkedList();
                RooCmdArg opt1 = RooFit::Save(); fitcmd->Add(&opt1);
                //RooCmdArg opt2 = RooFit::PrefitDataFraction(0.1); fitcmd->Add(&opt2);
                RooCmdArg opt3 = RooFit::Minimizer("Minuit","minimize"); fitcmd->Add(&opt3);//"migradimproved");//);
                RooCmdArg opt4 = RooFit::NumCPU(16, 0); fitcmd->Add(&opt4);
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
                drawText(Form("%i GeV < p_{T,#gamma#gamma} < %i GeV",pt,pt+1),0.6,0.68,1,16);
                drawText("|vz| < 30 cm",0.6,0.63,1,16);
                drawText("|#eta| < 1",0.6,0.58,1,16);
                drawText("p_{T,#gamma} > 1 GeV",0.6,0.53,1,16);
                drawText("#alpha < 0.6",0.6,0.48,1,16);

                const char * muname;
                const char * stdname;
                muname = "mean";
                stdname = "sigma1";
                RooRealVar* nsigVar = (RooRealVar*) res->floatParsFinal().find("nsig");
                RooRealVar* muVar = (RooRealVar*) res->floatParsFinal().find(muname);
                RooRealVar* stdVar = (RooRealVar*) res->floatParsFinal().find(stdname);
                float nsigyield = nsigVar->getVal();
                float nsigyielderr = nsigVar->getError(); 
                float mean = muVar->getVal();
                float meanerr = muVar->getError();
                float std = stdVar->getVal();
                float stderr = stdVar->getError(); 
                float chisq = 0;

                RooHist * hpull = plot->pullHist("data","model");
                double *ypull = hpull->GetY();
                int nFullBinsPull = 0;
                for(int i = 0; i < nBins; i++) {
                  if (ypull[i] == 0 || std::isnan(ypull[i])) continue;
                  chisq += TMath::Power(ypull[i],2);
                  nFullBinsPull++;
                }
                int numFitPar = res->floatParsFinal().getSize();
                int ndf = nFullBinsPull - numFitPar;
                drawText(Form("Fit #chi^2/NDF = %.2f",chisq/ndf),0.2,0.68,1,16);
                drawText(Form("#mu = %.4f#pm%.4f",mean,meanerr),0.2,0.63,1,16);
                drawText(Form("#sigma = %.4f#pm%.4f",std,stderr),0.2,0.58,1,16);
                drawText(Form("Yield: %.2f#pm%.2f",nsigyield,nsigyielderr),0.2,0.53,1,16);
                drawText(bkg->GetTitle(),0.2,0.48,1,16);
                c->cd();
                pplot->Draw();


                if (chisq/ndf < 6 && chisq/ndf > 0.05) {
                  cout << "success for " << pt << " " << a << " " << i << " " << j << " " << k << " " << l << endl;
                  fitpeaks[pt][a][i][j][k][l] = mean;
                  fitwidths[pt][a][i][j][k][l] = std;
                  fitchisq[pt][a][i][j][k][l] = chisq/ndf;
                  c->SaveAs("smearfits.pdf");
                  ws->Write();
                  delete h;
                  delete plot;
                  delete signal1;
                  delete model;
                  delete res;
                  delete pplot;
                  delete line;
                  delete poly2;
                  delete poly3;
                  delete ws;
                  break;
                }

                cout << "failure for " << pt << " " << a << " " << i << " " << j << " " << k << " " << l << endl;

                delete h;
                delete plot;
                delete signal1;
                delete model;
                delete res;
                delete pplot;
                delete line;
                delete poly2;
                delete poly3;
                delete ws;

              }
              wf->Close();
              delete wf;
            }
          }
        }
      }
    }
  }

  c->SaveAs("smearfits.pdf]");

  float bestratioscore = 100000000;
  int bestratio = 0;
  TH1D * hchisq = new TH1D("hchisq",";chisq;counts",1000,-5,200);
  for (int a = 0; a < 3; a++) {
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        for (int k = 0; k < 5; k++) {
          for (int l = 0; l < 5; l++) {
            float ratiodiff = 0;
            int count = 0;
            for (int pt = 0; pt < 3; pt++) {
              float chisq = fitchisq[pt][a][i][j][k][l];
              hchisq->Fill(chisq);
              if (chisq > 6 || chisq < 0.05) continue;
              ratiodiff += TMath::Power(datawidths[pt]/datapeaks[pt]-fitwidths[pt][a][i][j][k][l]/fitpeaks[pt][a][i][j][k][l],2);
              count++; 
            }
            if (count == 3) {
              ratiodiff = TMath::Sqrt(ratiodiff/(count));
              if (ratiodiff < bestratioscore) {
                bestratioscore = ratiodiff;
                bestratio = 1000*i + 100*j + 10*k + l;
              }
            }
          }}
      }
    }
  }
  cout << "bestratio: " << bestratio << " " << bestratioscore << endl;
  TFile * hfile = TFile::Open("smeartestchisq.root","RECREATE");
  hchisq->Write();
}

