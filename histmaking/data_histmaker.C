#include <iostream>
#include "ana.cxx"
#include "TreeSetting.h"

void data_histmaker(int runnum=48349)
{
  //gSystem->Load("../Headers/libMyDict.so");

  const char *dir = Form("/sphenix/tg/tg01/jets/samfred/run25/ppg11ana_hadded");
  string filename = Form("%s/run%i.root",dir,runnum);
  float ptcut = 1;
  bool doprob = false;

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
  //const char * wfilename = Form("test_mass_with_%i.root",runnum);
  TFile *wf = new TFile(wfilename,"recreate");
  
  // MBD trigger histos
  TH1D *hmass_mbd[nPtBins];
  TH1D *hmass_pi0_mbd[nPtBins];
  TH1D *hmass_eta_mbd[nPtBins];
  for(int ipt= 0; ipt < nPtBins; ipt++){
    hmass_mbd[ipt] = new TH1D(Form("hmass_pt%d_mbd",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
    hmass_pi0_mbd[ipt] = new TH1D(Form("hmass_pi0_pt%d_mbd",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
    hmass_eta_mbd[ipt] = new TH1D(Form("hmass_eta_pt%d_mbd",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  }
  TH1D * hmass_all_mbd = new TH1D(Form("hmass_mbd"),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
  TH1D * hmass_pi0_all_mbd = new TH1D(Form("hmass_pi0_mbd"),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
  TH1D * hmass_eta_all_mbd = new TH1D(Form("hmass_eta_mbd"),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  
  // photon3 trigger histos
  TH1D *hmass_photon3[nPtBins];
  TH1D *hmass_pi0_photon3[nPtBins];
  TH1D *hmass_eta_photon3[nPtBins];
  for(int ipt= 0; ipt < nPtBins; ipt++){
    hmass_photon3[ipt] = new TH1D(Form("hmass_pt%d_photon3",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
    hmass_pi0_photon3[ipt] = new TH1D(Form("hmass_pi0_pt%d_photon3",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
    hmass_eta_photon3[ipt] = new TH1D(Form("hmass_eta_pt%d_photon3",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  }
  TH1D * hmass_all_photon3 = new TH1D(Form("hmass_photon3"),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
  TH1D * hmass_pi0_all_photon3 = new TH1D(Form("hmass_pi0_photon3"),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
  TH1D * hmass_eta_all_photon3 = new TH1D(Form("hmass_eta_photon3"),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  
  // photon4 trigger histos
  TH1D *hmass_photon4[nPtBins];
  TH1D *hmass_pi0_photon4[nPtBins];
  TH1D *hmass_eta_photon4[nPtBins];
  for(int ipt= 0; ipt < nPtBins; ipt++){
    hmass_photon4[ipt] = new TH1D(Form("hmass_pt%d_photon4",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
    hmass_pi0_photon4[ipt] = new TH1D(Form("hmass_pi0_pt%d_photon4",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
    hmass_eta_photon4[ipt] = new TH1D(Form("hmass_eta_pt%d_photon4",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  }
  TH1D * hmass_all_photon4 = new TH1D(Form("hmass_photon4"),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
  TH1D * hmass_pi0_all_photon4 = new TH1D(Form("hmass_pi0_photon4"),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
  TH1D * hmass_eta_all_photon4 = new TH1D(Form("hmass_eta_photon4"),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  
  // photon5 trigger histos
  TH1D *hmass_photon5[nPtBins];
  TH1D *hmass_pi0_photon5[nPtBins];
  TH1D *hmass_eta_photon5[nPtBins];
  for(int ipt= 0; ipt < nPtBins; ipt++){
    hmass_photon5[ipt] = new TH1D(Form("hmass_pt%d_photon5",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
    hmass_pi0_photon5[ipt] = new TH1D(Form("hmass_pi0_pt%d_photon5",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
    hmass_eta_photon5[ipt] = new TH1D(Form("hmass_eta_pt%d_photon5",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  }
  TH1D * hmass_all_photon5 = new TH1D(Form("hmass_photon5"),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
  TH1D * hmass_pi0_all_photon5 = new TH1D(Form("hmass_pi0_photon5"),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
  TH1D * hmass_eta_all_photon5 = new TH1D(Form("hmass_eta_photon5"),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
    
  // Trigger efficiency
  TH1D * hmbd = new TH1D("hmbd",";p_{T}^{max};counts",1000,0,100);
  TH1D * hphoton3 = new TH1D("hphoton3",";p_{T}^{max};counts",1000,0,100);
  TH1D * hphoton4 = new TH1D("hphoton4",";p_{T}^{max};counts",1000,0,100);
  TH1D * hphoton5 = new TH1D("hphoton5",";p_{T}^{max};counts",1000,0,100);

  TH2D * hetaphi = new TH2D(Form("hetaphi"),";#eta;#phi",96,-1.1,1.1,256,-M_PI,M_PI);
  TH2D * hpi0_etaphi = new TH2D(Form("hpi0_etaphi"),";#eta;#phi",96,-1.1,1.1,256,-M_PI,M_PI);
  TH1D * hpt = new TH1D(Form("hpt"),";p_{T,#gamma};counts",20,0,20);

  Long64_t nentries = t->GetEntriesFast();
  Long64_t livecount=0;
  Long64_t scaledcount=0;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    t->GetEntry(jentry);
    if (fabs(vz) > 30) continue;
    bool foundmaxcluster = false;
    for(int icl=0;icl<nClusters_mother;icl++){
      auto pho = *(TLorentzVector*) photon_4mom_mother->At(icl);
      if(pho.E() > ptcut) foundmaxcluster = true;
    }
    if(!foundmaxcluster) continue;

    if(jentry % 1000==0) std::cout << "entry " << jentry << "/" << nentries << " (" << (float)jentry/nentries*100. << "%)" << std::endl;

    float maxpt = 0;
    for (int ic = 0; ic < nClusters; ic++) {
      auto pho = (TLorentzVector*) photon_4mom->At(ic);
      if (pho->E() > 1 && ScaledTriggerBit[10]/* && photon_prob[ic] > 0.05*/) hetaphi->Fill(pho->Eta(),pho->Phi());
      if (pho->Eta() < 1 && pho->Eta() > -1 && ScaledTriggerBit[10]/*&& photon_prob[ic] > 0.05*/) hpt->Fill(pho->Pt());
      if (pho->Pt() > maxpt) maxpt = pho->Pt();
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

      if(eta1 > etamax || eta1 < etamin) continue;
      if(eta2 > etamax || eta2 < etamin) continue;
      if(pt1 < ptcut || pt2 < ptcut) continue;

      if(diphoton_energyimbal[ip] > 0.6) continue;

      if (doprob) {
        if(photon_prob[idx_photon1[ip]] < 0.05 || photon_prob[idx_photon2[ip]] < 0.05) continue;
      }
      
      if (pho->M() < 0.16 && pho->M() > 0.11) hpi0_etaphi->Fill(pho->Eta(),pho->Phi());
      int bin = (int)pt;
      if (ScaledTriggerBit[10]) {
        if(bin >= 0 && bin < nPtBins){
          hmass_mbd[bin]->Fill(mass);
          hmass_pi0_mbd[bin]->Fill(mass);
          hmass_eta_mbd[bin]->Fill(mass);
        }
        hmass_all_mbd->Fill(mass);
        hmass_pi0_all_mbd->Fill(mass);
        hmass_eta_all_mbd->Fill(mass);
        hmbd->Fill(maxpt);
        if (LiveTriggerBit[25] || LiveTriggerBit[36]) {
          hphoton3->Fill(maxpt);
        }
        if (LiveTriggerBit[26] || LiveTriggerBit[37]) {
          hphoton4->Fill(maxpt);
        }
        if (LiveTriggerBit[27] || LiveTriggerBit[38]) {
          hphoton5->Fill(maxpt);
        }
      }
      if (ScaledTriggerBit[25] || ScaledTriggerBit[36]) {
        if(bin >= 0 && bin < nPtBins){
          hmass_photon3[bin]->Fill(mass);
          hmass_pi0_photon3[bin]->Fill(mass);
          hmass_eta_photon3[bin]->Fill(mass);
        }
        hmass_all_photon3->Fill(mass);
        hmass_pi0_all_photon3->Fill(mass);
        hmass_eta_all_photon3->Fill(mass);
      }
      if (ScaledTriggerBit[26] || ScaledTriggerBit[37]) {
        if(bin >= 0 && bin < nPtBins){
          hmass_photon4[bin]->Fill(mass);
          hmass_pi0_photon4[bin]->Fill(mass);
          hmass_eta_photon4[bin]->Fill(mass);
        }
        hmass_all_photon4->Fill(mass);
        hmass_pi0_all_photon4->Fill(mass);
        hmass_eta_all_photon4->Fill(mass);
      }
      if (ScaledTriggerBit[27] || ScaledTriggerBit[38]) {
        if(bin >= 0 && bin < nPtBins){
          hmass_photon5[bin]->Fill(mass);
          hmass_pi0_photon5[bin]->Fill(mass);
          hmass_eta_photon5[bin]->Fill(mass);
        }
        hmass_all_photon5->Fill(mass);
        hmass_pi0_all_photon5->Fill(mass);
        hmass_eta_all_photon5->Fill(mass);
      }
    }
  }

  wf->cd();
  for(int ip=0; ip < nPtBins; ip++){
    hmass_mbd[ip]->Write();
    hmass_pi0_mbd[ip]->Write();
    hmass_eta_mbd[ip]->Write();
    hmass_photon3[ip]->Write();
    hmass_pi0_photon3[ip]->Write();
    hmass_eta_photon3[ip]->Write();
    hmass_photon4[ip]->Write();
    hmass_pi0_photon4[ip]->Write();
    hmass_eta_photon4[ip]->Write();
    hmass_photon5[ip]->Write();
    hmass_pi0_photon5[ip]->Write();
    hmass_eta_photon5[ip]->Write();
  }
  hmass_all_mbd->Write();
  hmass_pi0_all_mbd->Write();
  hmass_eta_all_mbd->Write();
  hmass_all_photon3->Write();
  hmass_pi0_all_photon3->Write();
  hmass_eta_all_photon3->Write();
  hmass_all_photon4->Write();
  hmass_pi0_all_photon4->Write();
  hmass_eta_all_photon4->Write();
  hmass_all_photon5->Write();
  hmass_pi0_all_photon5->Write();
  hmass_eta_all_photon5->Write();

  hmbd->Write();
  hphoton3->Write();
  hphoton4->Write();
  hphoton5->Write();

  hetaphi->Write();
  hpi0_etaphi->Write();
  hpt->Write();
}
