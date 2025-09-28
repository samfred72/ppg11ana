using namespace RooFit;

void draw_yield(const char * particle = "pi0",const char * type = "data", bool smear = 0) {
  const char * formattedparticle;
  if (strcmp(particle,"pi0") == 0) formattedparticle="#pi^{0}";
  else formattedparticle="#eta";
  TH1D * pi0yields = new TH1D(Form("%syields",particle),Form(";p_{T,%s};Yield from fit",formattedparticle),20,0,20);
  const char * smearstr = "";
  if (smear) smearstr = "_smear";
  for (int i = 2; i <= 7; i++) {
    TFile * inf;
    if (strcmp(type,"data") == 0) inf = TFile::Open(Form("/sphenix/user/samfred/run25/ppg11/histmaking/sufficient_without_prob_1GeV/workspace_fits_%s_%s_pt%i_mbd.root",type,particle,i),"READ");
    else if (strcmp(type,"MC") == 0) inf = TFile::Open(Form("/sphenix/user/samfred/run25/ppg11/histmaking/sufficient_without_prob_1GeV%s/workspace_fits_%s_%s_pt%i_MB.root",smearstr,type,particle,i),"READ");
    if (!inf) continue;
    RooWorkspace * ws = (RooWorkspace*)inf->Get("workspace");
    RooRealVar * variable = ws->var(Form("nsig"));
    cout << variable << endl;
    if (variable) {
      
      float yield = variable->getVal();
      float yielderr = variable->getError();
      cout << i << " " << yield << " +/- " << yielderr << endl; 
      pi0yields->SetBinContent(i+1,yield);
      pi0yields->SetBinError(i+1,yielderr);
    }
  }
  gStyle->SetOptStat(0);
  pi0yields->SetMarkerStyle(20);
  pi0yields->SetMarkerColor(kRed);
  pi0yields->Draw();
  cout << "Writing file: " << Form("/sphenix/user/samfred/run25/ppg11/histmaking/yields_without_prob_1GeV%s/%s_%syields.root",smearstr,type,particle) << endl;
  TFile * wf = new TFile(Form("/sphenix/user/samfred/run25/ppg11/histmaking/yields_without_prob_1GeV%s/%s_%syields.root",smearstr,type,particle),"RECREATE");
  pi0yields->Write();
}
