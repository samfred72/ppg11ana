#include <iostream>
#include "ana.cxx"
#include "TreeSetting.h"

void truthhistmaker(int section = 0, const char * type = "MB")
{
  //gSystem->Load("../Headers/libMyDict.so");

  int runnum = 28; 
  const char *dir = Form("/sphenix/tg/tg01/jets/samfred/run25/pythia_%s_hadded",type);
  string filename = (section == -1) ? Form("%s/run%i_%s.root",dir,runnum,type) : Form("%s/run%i_%i.root",dir,runnum,section);
  float ptcut = 1;
  bool doprob = false;
  bool dosmear = true;
  // Scales for MC smearing
  float pfrac = 0.0004;
  float efrac = 0.04;
  float econst = 0.18;
  float escale = 1.026;

  TFile *f = new TFile(Form("%s",filename.c_str()),"read");
  TTree *t = (TTree*) f->Get("towerntup");

  treesetup(t,1);
   
  const double masslow = 0.0;
  const double masshigh = 1;

  ana anaclone;
  const int nPtBins = anaclone.nPtBins;
  const int nAlphaBins = anaclone.nAlphaBins;
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
  TH1D * ptss = new TH1D("ptss",";pt;counts",256,0,10);
  TH1D * ptds = new TH1D("ptds",";pt;counts",256,0,10);
  TH1D * ptsu = new TH1D("ptsu",";pt;counts",256,0,10);
  TH1D * ptdu = new TH1D("ptdu",";pt;counts",256,0,10);
  TH1D * hadist = new TH1D("hadist",";alpha;counts",100,0,1);
  TH1D * halpha[nPtBins][nAlphaBins];
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < nAlphaBins; j++) {
      halpha[i][j] = new TH1D(Form("halpha%i_%i",i,j),";mass;counts",100,0,0.3);
    }
  }
  TRandom rndm;
  TH1D *hmaxpt = new TH1D("hmaxpt",";max p_{T,#pi^0}^{truth};weighted counts;",64,0,64);
  Long64_t nentries = t->GetEntriesFast();
  Long64_t livecount=0;
  Long64_t scaledcount=0;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    t->GetEntry(jentry);

    if (fabs(vz) > 30) continue;
    float maxpt = 0;
    bool ispi0 = false;
    bool isineta = false;
    bool haspt = false;
    for (int it = 0; it < truthpar_n; it++) {
      auto pho = *(TLorentzVector*) truth_diphoton_4mom->At(it);
      if (truth_id[it] == 111 && pho.Pt() > maxpt) maxpt = pho.Pt();
      if (truth_id[it] == 111) ispi0 = true;
      if (truth_id[it] == 111 && abs(pho.Eta()) < 1) isineta = true; 
      if (truth_id[it] == 111 && pho.Pt() > 1) haspt = true;
    }
    hmaxpt->Fill(maxpt);  

    //if (!ispi0 || !isineta || !haspt) continue;
    //bool foundmaxcluster = false;
    //for(int icl=0;icl<nClusters_mother;icl++){
    //  auto pho = *(TLorentzVector*) photon_4mom_mother->At(icl);
    //  if(pho.E() > minClusterMaxE) foundmaxcluster = true;
    //}
    //if(!foundmaxcluster) continue;

    if(jentry % 1000==0) std::cout << "entry " << jentry << "/" << nentries << " (" << (float)jentry/nentries*100. << "%)" << std::endl;

    for(int ip=0; ip<nPairs; ip++){
      auto pho = (TLorentzVector*) diphoton_4mom->At(ip);
      auto pho1 = (TLorentzVector*) photon_4mom->At(idx_photon1[ip]);
      auto pho2 = (TLorentzVector*) photon_4mom->At(idx_photon2[ip]);

      float pt = pho->Pt();
      float e1 = pho1->E();
      float e2 = pho2->E();
      float pt1 = pho1->Pt();
      float pt2 = pho2->Pt();

      float etamin = anaclone.GetShiftedEta(vz,-1);
      float etamax = anaclone.GetShiftedEta(vz,1);
      
      float eta1 = pho1->Eta();
      float eta2 = pho2->Eta();
      float phi1 = pho1->Phi();
      float phi2 = pho2->Phi();
      
      float eta1fl = rndm.Gaus(eta1,pfrac/(e1)); 
      float phi1fl = rndm.Gaus(phi1,pfrac/(e1)); 
      float eta2fl = rndm.Gaus(eta2,pfrac/(e2)); 
      float phi2fl = rndm.Gaus(phi2,pfrac/(e2)); 

      float e1fl = rndm.Gaus(e1*escale,econst+efrac/TMath::Sqrt(e1));
      float e2fl = rndm.Gaus(e2*escale,econst+efrac/TMath::Sqrt(e2));
      float pt1fl = e1fl/cosh(eta1fl);
      float pt2fl = e2fl/cosh(eta2fl);

      bool fill = true;
      if (eta1 > etamax || eta1 < etamin) fill = false;
      if (eta2 > etamax || eta2 < etamin) fill = false;
      if (pho->M() < 0.1 || pho->M() > .2) fill = false;
      if (fill) {
        ptsu->Fill(pt1);
        ptsu->Fill(pt2);
        ptdu->Fill(pt);
      }


      if (dosmear) {
        if (eta1fl > etamax || eta1fl < etamin) continue;
        if (eta2fl > etamax || eta2fl < etamin) continue;
      }
      else {
        if (eta1 > etamax || eta1 < etamin) continue;
        if (eta2 > etamax || eta2 < etamin) continue;
      }
      TLorentzVector pho1fl;
      TLorentzVector pho2fl;
      pho1fl.SetPtEtaPhiE(pt1fl,eta1fl,phi1fl,e1fl);
      pho2fl.SetPtEtaPhiE(pt2fl,eta2fl,phi2fl,e2fl);
      TLorentzVector phofl = pho1fl + pho2fl;
      

      if (phofl.M() > 0.1 && phofl.M() < 0.2) {
        ptss->Fill(pt1fl);
        ptss->Fill(pt2fl);
        ptds->Fill(phofl.Pt());
      }
      if (dosmear) {
        if(pt1fl < ptcut || pt2fl < ptcut) continue;
      }
      else {
        if(pt1 < ptcut || pt2 < ptcut) continue;
      }
      if (doprob) {
        if(photon_prob[idx_photon1[ip]] < 0.05 || photon_prob[idx_photon2[ip]] < 0.05) continue;
      }
      float energyimbal = (dosmear ? abs(pho1fl.E() - pho2fl.E())/(pho1fl.E() + pho2fl.E()) : diphoton_energyimbal[ip]);
     
      float mass = (dosmear ? phofl.M() : pho->M()); 
      
      // alpha hists
      int bin = anaclone.findBin(pt);
      if (bin < 0) continue;
      int pbin = (bin < 4 ? 0 : 1);
      int abin = anaclone.findAlphaBin(energyimbal);
      hadist->Fill(energyimbal);
      if (abs(eta1) < 0.3 && abs(eta2) < 0.3) {
        halpha[pbin][abin]->Fill(mass);
      }
      if(energyimbal > 0.6) continue;

      
      if (dosmear) pt = phofl.Pt();
      
      hmass[bin]->Fill(mass);
      hmass_pi0[bin]->Fill(mass);
      hmass_eta[bin]->Fill(mass);
    }
  }

  wf->cd();
  for(int ip=0; ip < nPtBins; ip++){
    hmass[ip]->Write();
    hmass_pi0[ip]->Write();
    hmass_eta[ip]->Write();
  }
  hmaxpt->Write();
  ptss->Write();
  ptds->Write();
  ptsu->Write();
  ptdu->Write();
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 10; j++) {
      halpha[i][j]->Write();
    }
  }
  hadist->Write();
}
