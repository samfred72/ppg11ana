#include <iostream>
#include "ana.cxx"
#include "TreeSetting.h"

void data_histmaker(int runnum=48446)
{
  //gSystem->Load("../Headers/libMyDict.so");

  const char *dir = Form("/sphenix/tg/tg01/jets/samfred/run25/ppg11ana_hadded");
  string filename = Form("%s/run%i.root",dir,runnum);
  float minClusterMaxE = 0.5;

  TFile *f = new TFile(Form("%s",filename.c_str()),"read");
  TTree *t = (TTree*) f->Get("towerntup");

  treesetup(t,0);
   
  const double masslow = 0.0;
  const double masshigh = 1;

  ana anaclone;
  const int nPtBins = anaclone.nPtBins;
  double pi0mass = anaclone.pi0mass;
  double etamass = anaclone.etamass;

  const char * wfilename = Form("data_hists/mass_%i.root",runnum);
  TFile *wf = new TFile(wfilename,"recreate");
  TH1D *hmass[nPtBins];
  TH1D *hmass_pi0[nPtBins];
  TH1D *hmass_eta[nPtBins];
  for(int ipt= 0; ipt < nPtBins; ipt++){
    hmass[ipt] = new TH1D(Form("hmass_pt%d",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
    hmass_pi0[ipt] = new TH1D(Form("hmass_pi0_pt%d",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
    hmass_eta[ipt] = new TH1D(Form("hmass_eta_pt%d",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  }
  TH1D * hmass_all = new TH1D(Form("hmass"),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
  TH1D * hmass_pi0_all = new TH1D(Form("hmass_pi0"),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
  TH1D * hmass_eta_all = new TH1D(Form("hmass_eta"),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  TH2D * hetaphi = new TH2D(Form("hetaphi"),";#eta;#phi",96,-1.1,1.1,256,-M_PI,M_PI);
  TH2D * hpi0_etaphi = new TH2D(Form("hpi0_etaphi"),";#eta;#phi",96,-1.1,1.1,256,-M_PI,M_PI);
  TH1D * hpt = new TH1D(Form("hpt"),";p_{T,#gamma};counts",20,0,20);

  Long64_t nentries = t->GetEntriesFast();
  Long64_t livecount=0;
  Long64_t scaledcount=0;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    t->GetEntry(jentry);
    if (!ScaledTriggerBit[10]) continue;
    if (fabs(vz) > 30) continue;
    bool foundmaxcluster = false;
    for(int icl=0;icl<nClusters_mother;icl++){
      auto pho = *(TLorentzVector*) photon_4mom_mother->At(icl);
      if(pho.E() > minClusterMaxE) foundmaxcluster = true;
    }
    if(!foundmaxcluster) continue;

    if(jentry % 1000==0) std::cout << "entry " << jentry << "/" << nentries << " (" << (float)jentry/nentries*100. << "%)" << std::endl;

    for (int ic = 0; ic < nClusters; ic++) {
      auto pho = (TLorentzVector*) photon_4mom->At(ic);
      if (pho->E() > 0.5 && photon_prob[ic] > 0.05) hetaphi->Fill(pho->Eta(),pho->Phi());
      if (pho->Eta() < 1 && pho->Eta() > -1 && photon_prob[ic] > 0.05) hpt->Fill(pho->Pt());
    }

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
      
      if (pho->M() < 0.16 && pho->M() > 0.11) hpi0_etaphi->Fill(pho->Eta(),pho->Phi());
      int bin = (int)pt;
      if(bin >= 0 && bin < nPtBins){
        hmass[bin]->Fill(mass);
        hmass_pi0[bin]->Fill(mass);
        hmass_eta[bin]->Fill(mass);
      }
      hmass_all->Fill(mass);
      hmass_pi0_all->Fill(mass);
      hmass_eta_all->Fill(mass);
    }
  }

  wf->cd();
  for(int ip=0; ip < nPtBins; ip++){
    hmass[ip]->Write();
    hmass_pi0[ip]->Write();
    hmass_eta[ip]->Write();
  }
  hmass_all->Write();
  hmass_pi0_all->Write();
  hmass_eta_all->Write();
  hetaphi->Write();
  hpi0_etaphi->Write();
  hpt->Write();
}
