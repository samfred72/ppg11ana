#include <iostream>
#include "../Headers/ana.cxx"
#include "../Headers/TreeSetting.h"
#include <TEfficiency.h>

void MakeTruthHist(std::string particletype="pi0", std::string sampletype="MB", bool doReweight=false)
{
  int mbindex=3;
  gSystem->Load("../Headers/libMyDict.so");

  std::string trigstring;
  int triggerbit;
  int minClusterMaxE;

  ana anaclone;
  const int nbins = anaclone.nPtBins;
  int nbinsfine = anaclone.nPtBinsFine;
  const double *ptbin = anaclone.ptBins;//[nbins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  const double *ptbinfine = anaclone.ptBinsFine;
  double MAXETACUT = anaclone.MAXETACUT;

  const char *dir = Form("/sphenix/tg/tg01/jets/samfred/run25/pythia_%s_hadded/",sampletype.c_str());
  std::string filename;
  std::string rwstr = "noreweight";
  double masslow, masshigh, mass0;
  int pdgid=-1;
  
  TFile* fvzreweight;
  TH1D* hvzratio;
  TFile* fitreweight;
  TH1D* hreweight;
  if(doReweight){
    //fvzreweight = new TFile(Form("vz_reweight_%s.root",particletype.c_str()),"read");
    //hvzratio = (TH1D*) fvzreweight->Get("hratio");
    fitreweight = new TFile(Form("reweight_%s_%s.root",particletype.c_str(),sampletype.c_str()),"read");
    hreweight = (TH1D*) fitreweight->Get("hratio");
    rwstr = "reweight";
  }

  TFile *ftrigeff = new TFile("triggerefficiency.root","read");
  TF1* fitfunc = (TF1*) ftrigeff->Get("fitFunc");
  
  if(sampletype == "MB"){
    filename = Form("%s/run28_%d.root",dir,mbindex);
  }
  else if(sampletype== "Jet10"){
    filename = Form("%s/hadded_Run21_Pythia_Jet10.root",dir);
  }
  else if(sampletype == "Jet20"){
    filename = Form("%s/hadded_Run21_Pythia_Jet20.root",dir);
  }
  else if(sampletype == "Jet30"){
    filename = Form("%s/hadded_Run21_Pythia_Jet30.root",dir);
  }

  if(particletype == "pi0"){
    pdgid = 111;
  }
  else if(particletype == "eta"){
    pdgid = 221;
  }

  TFile *f = new TFile(Form("%s",filename.c_str()),"read");
  TTree *t = (TTree*) f->Get("towerntup");

  treesetup(t,1);
  const double vzcut = anaclone.vzcut;

  TFile *wf = new TFile(Form("output/hist_truth_%s_%s_%s_%d.root",rwstr.c_str(),particletype.c_str(), sampletype.c_str(),mbindex),"recreate");

  TH1D *heff_sp_pt_den = new TH1D("heff_sp_pt_den",";p_{T}^{truth} [GeV/c];",nbins,ptbin);
  TH1D *heff_sp_pt_den_fine = new TH1D("heff_sp_pt_den_fine",";p_{T}^{truth} [GeV/c];",nbinsfine,ptbinfine);
  TH1D *hmass= new TH1D("hmass",";mass;",50,0,1);
  TH1D *hmass_truth= new TH1D("hmass_truth",";mass;",50,0,1);
  
  Long64_t nentries = t->GetEntriesFast();
  Long64_t livecount=0;
  Long64_t scaledcount=0;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    t->GetEntry(jentry);
    float tretamin = anaclone.GetShiftedEta(truth_vz,-MAXETACUT);
    float tretamax = anaclone.GetShiftedEta(truth_vz,MAXETACUT);
    
    float truthpt = -1;
    float truthy = -999;
    float truthphi = -999;
    float trutheta = -999;
    float trutheta1 = -999;
    float trutheta2 = -999;
    float truthphi1 = -999;
    float truthphi2 = -999;
    float truthe = -1;

    int maxtruthidx = -1;
    bool foundmaxpt = false;
    double vzreweightfactor = 1.;//(doReweight) ? hvzratio->GetBinContent(hvzratio->FindFixBin(truth_vz)) : 1;
    double reweightfactor = 1;
    int ntr = 0;
    if(jentry % 5000==0) std::cout << "entry " << jentry << "/" << nentries << " (" << (float)jentry/nentries*100. << "%)" << std::endl;
    for(int itruth=0;itruth < truthpar_n; itruth++){
      if(fabs(truth_vz)>vzcut)continue;
      if(truth_id[itruth] != pdgid) continue;
      auto trdipho = *(TLorentzVector*) truth_diphoton_4mom->At(itruth);
      auto trpho1 = *(TLorentzVector*) truth_photon1_4mom->At(itruth);
      auto trpho2 = *(TLorentzVector*) truth_photon2_4mom->At(itruth);
      float _truthpt = trdipho.Pt();
      float _truthy = trdipho.Rapidity();
      if(_truthpt<1) continue;
      if(_truthy < tretamin || _truthy > tretamax ) continue;
      if(particletype=="pi0" && (!truth_found_decay1[itruth] || !truth_found_decay2[itruth])) continue;
      //if(!truth_found_decay1[itruth] || !truth_found_decay2[itruth]) continue;
      if(_truthpt > truthpt){
        truthpt = _truthpt;
        truthy = _truthy;
        trutheta = trdipho.Eta();
        trutheta1 = trpho1.Eta();
        trutheta2 = trpho2.Eta();
        truthphi = trdipho.Phi();
        truthphi1 = trpho1.Phi();
        truthphi2 = trpho2.Phi();
        truthe = trdipho.E();
        maxtruthidx = itruth;
        foundmaxpt = true; 
      }
      reweightfactor =  (doReweight) ? hreweight->GetBinContent(hreweight->FindFixBin(_truthpt)) * vzreweightfactor : 1;
      heff_sp_pt_den->Fill(_truthpt,reweightfactor);
      heff_sp_pt_den_fine->Fill(_truthpt,reweightfactor);
      ntr++;
    } 
    
    if(fabs(vz)>30) continue;
    for(int ip=0; ip<nPairs; ip++){
      auto pho = (TLorentzVector*) diphoton_4mom->At(ip);
      auto pho1 = (TLorentzVector*) photon_4mom->At(idx_photon1[ip]);
      auto pho2 = (TLorentzVector*) photon_4mom->At(idx_photon2[ip]);

      float pt1 = pho1->Pt();
      float pt2 = pho2->Pt();

      float en = pho->E();
      float pt = pho->Pt();
      if(pt<2 || pt>3) continue;
      float y = pho->Rapidity();
      float eta = pho->Eta();
      float phi = pho->Phi();
      float mass = pho->M();
      float etamin = anaclone.GetShiftedEta(vz,-MAXETACUT);
      float etamax = anaclone.GetShiftedEta(vz,MAXETACUT);
      
      float eta1 = pho1->Eta();
      float eta2 = pho2->Eta();
      float phi1 = pho1->Phi();
      float phi2 = pho2->Phi();

      if(eta1>etamax || eta1 < etamin) continue;
      if(eta2>etamax || eta2 < etamin) continue;
      if(pt1<0.5 || pt2<0.5) continue;
      if(diphoton_energyimbal[ip]>0.6) continue;

      if(photon_prob[idx_photon1[ip]] < 0.05) continue;
      if(photon_prob[idx_photon2[ip]] < 0.05) continue;
      hmass->Fill(mass);
      if(ntr>0) hmass_truth->Fill(mass);
    }
  }

  wf->cd();
  heff_sp_pt_den->Write();
  heff_sp_pt_den_fine->Write();
  hmass->Write();
  hmass_truth->Write();
}
