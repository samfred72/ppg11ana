#include <iostream>
#include "ana.cxx"
#include "TreeSetting.h"

void histmaker(int section = 0, const char * type = "MB")
{
  gSystem->Load("../Headers/libMyDict.so");

  int runnum = 28;
  const char *dir = Form("/sphenix/tg/tg01/jets/samfred/run25/pythia_%s_hadded",type);
  string filename = (section == -1) ? Form("%s/run%i_%s.root",dir,runnum,type) : Form("%s/run%i_%i.root",dir,runnum,section);
  float minClusterMaxE = 0.5;

  TFile *f = new TFile(Form("%s",filename.c_str()),"read");
  TTree *t = (TTree*) f->Get("towerntup");

  treesetup(t,1);
   
  const double masslow = 0.0;
  const double masshigh = 1;

  ana anaclone;
  const int nPtBins = anaclone.nPtBins;
  double pi0mass = anaclone.pi0mass;
  double etamass = anaclone.etamass;

  const char * wfilename = (section == -1) ? Form("hists/mass_%s.root",type) : Form("hists/mass_%s_%i.root",type,section);
  TFile *wf = new TFile(wfilename,"recreate");
  TH1D *hmass[nPtBins];
  TH1D *hmass_pi0[nPtBins];
  TH1D *hmass_eta[nPtBins];
  for(int ipt= 0; ipt < nPtBins; ipt++){
    hmass[ipt] = new TH1D(Form("hmass_pt%d",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
    hmass_pi0[ipt] = new TH1D(Form("hmass_pi0_pt%d",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
    hmass_eta[ipt] = new TH1D(Form("hmass_eta_pt%d",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  }
  TH1D *hmaxpt = new TH1D("hmaxpt",";max p_{T,#pi^0}^{truth};weighted counts;",64,0,64);
  Long64_t nentries = t->GetEntriesFast();
  Long64_t livecount=0;
  Long64_t scaledcount=0;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    t->GetEntry(jentry);

    if (fabs(vz) > 30) continue;
    float maxpt = 0;
    for (int it = 0; it < truthpar_n; it++) {
      auto pho = *(TLorentzVector*) truth_diphoton_4mom->At(it);
      if (truth_id[it] == 111 && pho.Pt() > maxpt) maxpt = pho.Pt();
    }
    hmaxpt->Fill(maxpt);  

    bool foundmaxcluster = false;
    for(int icl=0;icl<nClusters_mother;icl++){
      auto pho = *(TLorentzVector*) photon_4mom_mother->At(icl);
      if(pho.E() > minClusterMaxE) foundmaxcluster = true;
    }
    if(!foundmaxcluster) continue;

    if(jentry % 1000==0) std::cout << "entry " << jentry << "/" << nentries << " (" << (float)jentry/nentries*100. << "%)" << std::endl;

    for(int ip=0; ip<nPairs; ip++){
      auto pho = (TLorentzVector*) diphoton_4mom->At(ip);
      auto pho1 = (TLorentzVector*) photon_4mom->At(idx_photon1[ip]);
      auto pho2 = (TLorentzVector*) photon_4mom->At(idx_photon2[ip]);

      float pt1 = pho1->Pt();
      float pt2 = pho2->Pt();

      float pt = pho->Pt();
      float mass = pho->M();
      float etamin = anaclone.GetShiftedEta(vz,-1);
      float etamax = anaclone.GetShiftedEta(vz,1);
      
      float eta1 = pho1->Eta();
      float eta2 = pho2->Eta();
      float phi1 = pho1->Phi();
      float phi2 = pho2->Phi();

      if(fabs(eta1)>etamax || fabs(eta1) < etamin) continue;
      if(fabs(eta2)>etamax || fabs(eta2) < etamin) continue;
      if(pt1 < 0.5 || pt2 < 0.5) continue;

      if(diphoton_energyimbal[ip] > 0.6) continue;

      if(photon_prob[idx_photon1[ip]] < 0.05 || photon_prob[idx_photon2[ip]] < 0.05) continue;
      
      int bin = (int)pt;
      if(bin >= 0 && bin < nPtBins){
        hmass[bin]->Fill(mass);
        hmass_pi0[bin]->Fill(mass);
        hmass_eta[bin]->Fill(mass);
      }
    }
  }

  wf->cd();
  for(int ip=0; ip < nPtBins; ip++){
    hmass[ip]->Write();
    hmass_pi0[ip]->Write();
    hmass_eta[ip]->Write();
  }
  hmaxpt->Write();
}
