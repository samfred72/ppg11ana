#include <iostream>
#include "ana.cxx"
#include "TreeSetting.h"

void toyMC(int section = 0, const char * type = "MB")
{
  gSystem->Load("../Headers/libMyDict.so");

  int runnum = 28; 
  const char *dir = Form("/sphenix/tg/tg01/jets/samfred/run25/pythia_%s_hadded",type);
  string filename = (section == -1) ? Form("%s/run%i_%s.root",dir,runnum,type) : Form("%s/run%i_%i.root",dir,runnum,section);

  TFile *f = new TFile(Form("%s",filename.c_str()),"read");
  TTree *t = (TTree*) f->Get("towerntup");

  treesetup(t,1);
 
  float esmearfactor = 0.1;
  float tsmearfactor = 0.001;
  TH1D * hist[10][10];
  for (int i = 0; i < 10; i++ ) {
    for (int j = 0; j < 10; j++) {
      hist[i][j] = new TH1D(Form("hist%i%i",i,j),Form("Theta smear %f, E smear %f;Mass [GeV];counts",i*tsmearfactor,j*esmearfactor),1000,0,0.3);
    }
  }
  TH1D * dtheta = new TH1D("dtheta",";#Delta#theta;counts",100,-M_PI,M_PI);
  int nentries = t->GetEntries()/10;
  TRandom rand;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    t->GetEntry(jentry);
    if (jentry % 100000 == 0) cout << "Processing event: (" << jentry << "/" << nentries << "): " << (int)((float)jentry/(float)nentries*100) << "%" << endl;
    for (int it = 0; it < truth_photon1_4mom->GetEntries(); it++) {
      TLorentzVector pho1 = *(TLorentzVector*) truth_photon1_4mom->At(it);
      TLorentzVector pho2 = *(TLorentzVector*) truth_photon2_4mom->At(it);
      TLorentzVector dipho = pho1 + pho2;
      float eta1 = pho1.Eta();
      float eta2 = pho2.Eta();
      float phi1 = pho1.Phi();
      float phi2 = pho2.Phi();
      float pt1 = pho1.Pt();
      float pt2 = pho2.Pt();
      float e1 = pho1.E();
      float e2 = pho2.E();
      if (abs(eta1) > 1 || abs(eta2) > 1) continue;
      if (pt1 < 0.5 || pt2 < 0.5) continue;
      float theta1 = 2*TMath::ATan(TMath::Exp(-eta1));
      float theta2 = 2*TMath::ATan(TMath::Exp(-eta2));
      float theta = pho1.Angle(pho2.Vect());
      dtheta->Fill(theta);
      for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++ ) {
          float smeartheta1 = rand.Gaus(theta1,0.0001*i);
          float smeartheta2 = rand.Gaus(theta2,0.0001*i);
          //float smeareta1 = -TMath::Log(TMath::Tan(smeartheta1/2.0));
          //float smeareta2 = -TMath::Log(TMath::Tan(smeartheta2/2.0));
          float smeareta1 = rand.Gaus(eta1,0.0001*i);
          float smeareta2 = rand.Gaus(eta2,0.0001*i);
          float smearphi1 = rand.Gaus(phi1,0.0001*i);
          float smearphi2 = rand.Gaus(phi2,0.0001*i);
          float smeare1 = rand.Gaus(e1,esmearfactor*j*e1);
          float smeare2 = rand.Gaus(e2,esmearfactor*j*e2);
          float smearpt1 = smeare1/TMath::CosH(eta1);
          float smearpt2 = smeare2/TMath::CosH(eta2);
          TLorentzVector smearpho1;
          TLorentzVector smearpho2;
          smearpho1.SetPtEtaPhiE(smearpt1,smeareta1,smearphi1,smeare1);
          smearpho2.SetPtEtaPhiE(smearpt2,smeareta2,smearphi2,smeare2);
          //smearpho1.SetPtEtaPhiE(pt1,smeareta1,phi1,e1);
          //smearpho2.SetPtEtaPhiE(pt2,eta2,phi2,e2);
          //smearpho1.SetPtEtaPhiE(smearpt1,eta1,phi1,smeare1);
          //smearpho2.SetPtEtaPhiE(smearpt2,eta2,phi2,smeare2);
          TLorentzVector smeardipho = smearpho1+smearpho2;
          //float smearedtheta = rand.Gaus(theta,tsmearfactor*i);
          //float mass = TMath::Sqrt(2*smeare1*smeare2*(1-TMath::Cos(smearedtheta)));
          if (smeardipho.Pt() < 2 || smeardipho.Pt() > 3) continue;
          hist[i][j]->Fill(smeardipho.M());
          //hist[i][j]->Fill(mass);
        }
      }
    }
  }
  TFile * wf = TFile::Open("toymass.root","RECREATE");
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 10; j++) {
      hist[i][j]->Write();
    }
  }
  dtheta->Write();
}
