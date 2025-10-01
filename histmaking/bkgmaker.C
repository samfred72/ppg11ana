#include <iostream>
#include "ana.cxx"
#include "TreeSetting.h"

void bkgmaker(int runnum=48349)
{
  //gSystem->Load("../Headers/libMyDict.so");

  const char *dir = Form("/sphenix/tg/tg01/jets/samfred/run25/ppg11ana_hadded");
  string filename = Form("%s/run%i.root",dir,runnum);
  float ptcut = 0.5;
  bool doprob = true;

  TFile *f = new TFile(Form("%s",filename.c_str()),"read");
  TTree *t = (TTree*) f->Get("towerntup");

  treesetup(t,0);
   
  const double masslow = 0.0;
  const double masshigh = 1;

  ana anaclone;
  const int nPtBins = anaclone.nPtBins;
  double pi0mass = anaclone.pi0mass;
  double etamass = anaclone.etamass;

  const char * wfilename = Form("bkg_hists/bkg_%i.root",runnum);
  //const char * wfilename = Form("test_mass_with_%i.root",runnum);
  TFile *wf = new TFile(wfilename,"recreate");
  
  // MBD trigger histos
  TH1D *hbkg_mbd[nPtBins];
  TH1D *hbkg_pi0_mbd[nPtBins];
  TH1D *hbkg_eta_mbd[nPtBins];
  for(int ipt= 0; ipt < nPtBins; ipt++){
    hbkg_mbd[ipt] = new TH1D(Form("hbkg_pt%d_mbd",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
    hbkg_pi0_mbd[ipt] = new TH1D(Form("hbkg_pi0_pt%d_mbd",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
    hbkg_eta_mbd[ipt] = new TH1D(Form("hbkg_eta_pt%d_mbd",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  }
  TH1D * hbkg_all_mbd = new TH1D(Form("hbkg_mbd"),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
  TH1D * hbkg_pi0_all_mbd = new TH1D(Form("hbkg_pi0_mbd"),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
  TH1D * hbkg_eta_all_mbd = new TH1D(Form("hbkg_eta_mbd"),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  
  // photon3 trigger histos
  TH1D *hbkg_photon3[nPtBins];
  TH1D *hbkg_pi0_photon3[nPtBins];
  TH1D *hbkg_eta_photon3[nPtBins];
  for(int ipt= 0; ipt < nPtBins; ipt++){
    hbkg_photon3[ipt] = new TH1D(Form("hbkg_pt%d_photon3",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
    hbkg_pi0_photon3[ipt] = new TH1D(Form("hbkg_pi0_pt%d_photon3",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
    hbkg_eta_photon3[ipt] = new TH1D(Form("hbkg_eta_pt%d_photon3",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  }
  TH1D * hbkg_all_photon3 = new TH1D(Form("hbkg_photon3"),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
  TH1D * hbkg_pi0_all_photon3 = new TH1D(Form("hbkg_pi0_photon3"),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
  TH1D * hbkg_eta_all_photon3 = new TH1D(Form("hbkg_eta_photon3"),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  
  // photon4 trigger histos
  TH1D *hbkg_photon4[nPtBins];
  TH1D *hbkg_pi0_photon4[nPtBins];
  TH1D *hbkg_eta_photon4[nPtBins];
  for(int ipt= 0; ipt < nPtBins; ipt++){
    hbkg_photon4[ipt] = new TH1D(Form("hbkg_pt%d_photon4",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
    hbkg_pi0_photon4[ipt] = new TH1D(Form("hbkg_pi0_pt%d_photon4",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
    hbkg_eta_photon4[ipt] = new TH1D(Form("hbkg_eta_pt%d_photon4",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  }
  TH1D * hbkg_all_photon4 = new TH1D(Form("hbkg_photon4"),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
  TH1D * hbkg_pi0_all_photon4 = new TH1D(Form("hbkg_pi0_photon4"),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
  TH1D * hbkg_eta_all_photon4 = new TH1D(Form("hbkg_eta_photon4"),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  
  // photon5 trigger histos
  TH1D *hbkg_photon5[nPtBins];
  TH1D *hbkg_pi0_photon5[nPtBins];
  TH1D *hbkg_eta_photon5[nPtBins];
  for(int ipt= 0; ipt < nPtBins; ipt++){
    hbkg_photon5[ipt] = new TH1D(Form("hbkg_pt%d_photon5",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
    hbkg_pi0_photon5[ipt] = new TH1D(Form("hbkg_pi0_pt%d_photon5",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
    hbkg_eta_photon5[ipt] = new TH1D(Form("hbkg_eta_pt%d_photon5",ipt),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  }
  TH1D * hbkg_all_photon5 = new TH1D(Form("hbkg_photon5"),";m_{#gamma_{1}#gamma_{2}};",100,0,1);
  TH1D * hbkg_pi0_all_photon5 = new TH1D(Form("hbkg_pi0_photon5"),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
  TH1D * hbkg_eta_all_photon5 = new TH1D(Form("hbkg_eta_photon5"),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
  

  TH1D * hpt = new TH1D(Form("hpt"),";p_{T,#gamma};counts",20,0,20);

  bool hasmbdphoton = false;
  bool hasphoton3photon = false;
  bool hasphoton4photon = false;
  bool hasphoton5photon = false;
  bool hasnextmbdphoton = false;
  bool hasnextphoton3photon = false;
  bool hasnextphoton4photon = false;
  bool hasnextphoton5photon = false;
  TLorentzVector mbdpho; 
  TLorentzVector photon3pho; 
  TLorentzVector photon4pho; 
  TLorentzVector photon5pho; 
  TLorentzVector nextmbdpho; 
  TLorentzVector nextphoton3pho; 
  TLorentzVector nextphoton4pho; 
  TLorentzVector nextphoton5pho; 

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

    for (int ic = 0; ic < nClusters; ic++) {
      auto pho = (TLorentzVector*) photon_4mom->At(ic);
      TLorentzVector diphoton = *pho + mbdpho;
      float pt = diphoton.Pt();
      float mass = diphoton.M();

      float etamin = anaclone.GetShiftedEta(vz,-1);
      float etamax = anaclone.GetShiftedEta(vz,1);
      
      if (eta1 > etamax || eta1 < etamin) continue;
      if (eta2 > etamax || eta2 < etamin) continue;
      if (pho->Pt() < ptcut) continue;
      if (doprob) {
        if(photon_prob[ic] < 0.05) continue;
      }

      if (ScaledTriggerBit[10] && !hasnextmbdphoton) {
        nextmbdphoton = pho;
        hasnextmbdphoton = true;
      }
      if ((ScaledTriggerBit[25] || ScaledTriggerBit[36]) && !hasnextphoton3photon) {
        nextphoton3photon = pho;
        hasnextphoton3photon = true;
      }
      if ((ScaledTriggerBit[26] || ScaledTriggerBit[37]) && !hasnextphoton4photon) {
        nextphoton4photon = pho;
        hasnextphoton4photon = true;
      }
      if ((ScaledTriggerBit[27] || ScaledTriggerBit[38]) && !hasnextphoton5photon) {
        nextphoton5photon = pho;
        hasnextphoton5photon = true;
      }

      int bin = (int)pt;
      if (ScaledTriggerBit[10]) {
        if (!hasmbdphoton) {
          mbdphoton = pho; 
          continue;
        }
        if(bin >= 0 && bin < nPtBins){
          hbkg_mbd[bin]->Fill(mass);
          hbkg_pi0_mbd[bin]->Fill(mass);
          hbkg_eta_mbd[bin]->Fill(mass);
        }
        hbkg_all_mbd->Fill(mass);
        hbkg_pi0_all_mbd->Fill(mass);
        hbkg_eta_all_mbd->Fill(mass);
      }
      if (ScaledTriggerBit[25] || ScaledTriggerBit[36]) {
        if (!hasphoton3photon) {
          photon3photon = pho; 
          continue;
        }
        if(bin >= 0 && bin < nPtBins){
          hbkg_photon3[bin]->Fill(mass);
          hbkg_pi0_photon3[bin]->Fill(mass);
          hbkg_eta_photon3[bin]->Fill(mass);
        }
        hbkg_all_photon3->Fill(mass);
        hbkg_pi0_all_photon3->Fill(mass);
        hbkg_eta_all_photon3->Fill(mass);
      }
      if (ScaledTriggerBit[26] || ScaledTriggerBit[37]) {
        if (!hasphoton4photon) {
          photon4photon = pho; 
          continue;
        }
        if(bin >= 0 && bin < nPtBins){
          hbkg_photon4[bin]->Fill(mass);
          hbkg_pi0_photon4[bin]->Fill(mass);
          hbkg_eta_photon4[bin]->Fill(mass);
        }
        hbkg_all_photon4->Fill(mass);
        hbkg_pi0_all_photon4->Fill(mass);
        hbkg_eta_all_photon4->Fill(mass);
      }
      if (ScaledTriggerBit[27] || ScaledTriggerBit[38]) {
        if (!hasphoton5photon) {
          photon5photon = pho; 
          continue;
        }
        if(bin >= 0 && bin < nPtBins){
          hbkg_photon5[bin]->Fill(mass);
          hbkg_pi0_photon5[bin]->Fill(mass);
          hbkg_eta_photon5[bin]->Fill(mass);
        }
        hbkg_all_photon5->Fill(mass);
        hbkg_pi0_all_photon5->Fill(mass);
        hbkg_eta_all_photon5->Fill(mass);
      }

    }

    for(int ip=0; ip<nPairs; ip++){
      auto pho = (TLorentzVector*) diphoton_4mom->At(ip);
      auto pho1 = (TLorentzVector*) photon_4mom->At(idx_photon1[ip]);
      auto pho2 = (TLorentzVector*) photon_4mom->At(idx_photon2[ip]);

      float pt1 = pho1->Pt();
      float pt2 = pho2->Pt();

      float pt = pho->Pt();
      float mass = pho->M();
      
      float eta1 = pho1->Eta();
      float eta2 = pho2->Eta();
      float phi1 = pho1->Phi();
      float phi2 = pho2->Phi();

      if(fabs(eta1)>etamax || fabs(eta1) < etamin) continue;
      if(fabs(eta2)>etamax || fabs(eta2) < etamin) continue;
      if(pt1 < ptcut || pt2 < ptcut) continue;

      if(diphoton_energyimbal[ip] > 0.6) continue;

      
      if (pho->M() < 0.16 && pho->M() > 0.11) hpi0_etaphi->Fill(pho->Eta(),pho->Phi());
    }
  }

  wf->cd();
  for(int ip=0; ip < nPtBins; ip++){
    hbkg_mbd[ip]->Write();
    hbkg_pi0_mbd[ip]->Write();
    hbkg_eta_mbd[ip]->Write();
    hbkg_photon3[ip]->Write();
    hbkg_pi0_photon3[ip]->Write();
    hbkg_eta_photon3[ip]->Write();
    hbkg_photon4[ip]->Write();
    hbkg_pi0_photon4[ip]->Write();
    hbkg_eta_photon4[ip]->Write();
    hbkg_photon5[ip]->Write();
    hbkg_pi0_photon5[ip]->Write();
    hbkg_eta_photon5[ip]->Write();
  }
  hbkg_all_mbd->Write();
  hbkg_pi0_all_mbd->Write();
  hbkg_eta_all_mbd->Write();
  hbkg_all_photon3->Write();
  hbkg_pi0_all_photon3->Write();
  hbkg_eta_all_photon3->Write();
  hbkg_all_photon4->Write();
  hbkg_pi0_all_photon4->Write();
  hbkg_eta_all_photon4->Write();
  hbkg_all_photon5->Write();
  hbkg_pi0_all_photon5->Write();
  hbkg_eta_all_photon5->Write();

  hpt->Write();
}
