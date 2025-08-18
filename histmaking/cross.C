void cross(const char * particle = "pi0") {
  TFile * fyield = TFile::Open(Form("/sphenix/user/samfred/run25/ppg11/%syields.root",particle),"READ");
  TFile * faccept = TFile::Open(Form("/sphenix/user/jpark4/Public/ForSam/eff_output_noreweight_pythiaMB_%s.root",particle),"READ");
  const char * fitfilename;
  if (strcmp(particle,"pi0") == 0) fitfilename = "fit_pythiamb_cross_section_pi0.root";
  else fitfilename = "fit_pythiamb_cross_section_tsallis_eta.root";
  TFile * fave = TFile::Open(Form("/sphenix/user/jpark4/Public/ForSam/%s",fitfilename),"READ");
  TH1D * hyield = (TH1D*)fyield->Get(Form("%syields",particle));
  TH1D * haccept = (TH1D*)faccept->Get(Form("heff"));
  TF1 * ave = (TF1*)fave->Get("fit_restricted");
  TF1 * ave_x_pt = (TF1*)fave->Get("f_pTdNdpt");

  const char * formattedparticle;
  if (strcmp(particle,"pi0") == 0) formattedparticle="#pi^{0}";
  TFile * wf = new TFile(Form("cross_section_%s.root",particle),"RECREATE");
  TH1D * hcross = new TH1D("hcross",Form(";p_{T,%s};E#frac{d^{3}#sigma}{dp^{3}} [mb GeV^{-2} c^{3}",formattedparticle),20,0,20);
  float L = 100e6/25.2;
  float B = 0.98823;
  float dy = 2;
  for (int i = 1; i < hcross->GetNbinsX(); i++) {
    float binlow = hcross->GetXaxis()->GetBinLowEdge(i);
    float binup = hcross->GetXaxis()->GetBinUpEdge(i);
    float avept = ave_x_pt->Integral(binlow,binup)/ave->Integral(binlow,binup);
    float deltapt = binup-binlow;
    float yeild = hyield->GetBinContent(i);
    float accept = haccept->GetBinContent(i);
    if (avept == 0 || deltapt == 0 || deltay == 0 || L == 0 || B == 0 || accept == 0) continue; 
    float cross = yeild/avept/accept/deltapt;
    hcross->SetBinContent(i,cross);
    float yeilderr = hyield->GetBinError(i);
    hcross->SetBinError(i,yeilderr);
  }
  hcross->Scale(1.0/2/M_PI/L/B/dy);
  
  for (int i = 0; i < hcross->GetNbinsX(); i++) {
    cout << hcross->GetBinContent(i) << endl;
  }

  TCanvas * c = new TCanvas("c","",500,500);
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  hcross->SetMarkerStyle(20);
  hcross->SetMarkerColor(kRed);
  hcross->Draw("p");
}



