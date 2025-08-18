using namespace RooFit;

void draw_yield(const char * particle = "pi0") {
  TFile * inf = TFile::Open(Form("sufficient/workspace_fits_MB_%s.root",particle),"READ");
  RooWorkspace * ws = (RooWorkspace*)inf->Get("workspace");
  ws->Print();
  const char * formattedparticle;
  if (strcmp(particle,"pi0") == 0) formattedparticle="#pi^{0}";
  else formattedparticle="#eta";
  TH1D * pi0yields = new TH1D(Form("%syields",particle),Form(";p_{T,%s};Yield from fit",formattedparticle),20,0,20);
  for (int i = 0; i < 20; i++) {
    RooRealVar * variable = ws->var(Form("nsig%i",i));
    cout << variable << endl;
    if (variable) {
      
      float yield = variable->getVal();
      float yielderr = variable->getError();
      cout << yield << " +/- " << yielderr << endl; 
      pi0yields->SetBinContent(i+1,yield);
      pi0yields->SetBinError(i+1,yielderr);
    }
  }
  gStyle->SetOptStat(0);
  pi0yields->SetMarkerStyle(20);
  pi0yields->SetMarkerColor(kRed);
  pi0yields->Draw();
  TFile * wf = new TFile(Form("%syields.root",particle),"RECREATE");
  pi0yields->Write();
}
