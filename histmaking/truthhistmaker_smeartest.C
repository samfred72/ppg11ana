#include <iostream>
#include "ana.cxx"
#include "TreeSetting.h"

void truthhistmaker_smeartest(int section = 0, const char * type = "MB")
{
  //gSystem->Load("../Headers/libMyDict.so");

  ana anaclone;
  int runnum = 28; 
  const char *dir = Form("/sphenix/tg/tg01/jets/samfred/run25/pythia_%s_hadded",type);
  string filename = (section == -1) ? Form("%s/run%i_%s.root",dir,runnum,type) : Form("%s/run%i_%i.root",dir,runnum,section);
  int smears = 5;
  float ptcut = 1;
  bool doprob = false;
  bool dosmear = true;
  // Scales for MC smearing
  float pfrac = 0.001; float pfrac_high = anaclone.pfrac_high;  float pfrac_low = anaclone.pfrac_low;
  float efrac = 0.1; float efrac_high =   anaclone.efrac_high;  float efrac_low = anaclone.efrac_low;
  float econst = 0.1; float econst_high = anaclone.econst_high; float econst_low =anaclone.econst_low;
  float escale = 1.03; float escale_high =anaclone.escale_high; float escale_low =anaclone.escale_low;

  TFile *f = new TFile(Form("%s",filename.c_str()),"read");
  TTree *t = (TTree*) f->Get("towerntup");

  treesetup(t,1);
   
  const double masslow = 0.0;
  const double masshigh = 1;

  const int nPtBins = 3;//anaclone.nPtBins;
  const int nAlphaBins = 3;//anaclone.nAlphaBins;
  double pi0mass = anaclone.pi0mass;
  double etamass = anaclone.etamass;

  const char * wfilename = (section == -1) ? Form("hists/mass_%s.root",type) : Form("hists/mass_%s_%i.root",type,section);
  TFile *wf = new TFile(wfilename,"recreate");
  TH1D *hmass_pi0[nPtBins][nAlphaBins][smears][smears][smears][smears];
  TH1D *hmass_eta[nPtBins][nAlphaBins][smears][smears][smears][smears];
  for(int ipt= 0; ipt < nPtBins; ipt++) {
    for (int ia = 0; ia < nAlphaBins; ia++) {
      for (int i = 0; i < smears; i++) {
        for (int j = 0; j < smears; j++) {
          for (int k = 0; k < smears; k++) {
            for (int l = 0; l < smears; l++) {
              hmass_pi0[ipt][ia][i][j][k][l] = new TH1D(Form("hmass_pi0_pt%d_a%i_%i_%i_%i_%i",ipt,ia,i,j,k,l),";m_{#gamma_{1}#gamma_{2}};",100,0,0.3);
              hmass_eta[ipt][ia][i][j][k][l] = new TH1D(Form("hmass_eta_pt%d_a%i_%i_%i_%i_%i",ipt,ia,i,j,k,l),";m_{#gamma_{1}#gamma_{2}};",100,0.3,0.8);
            }
          }
        }
      }
    }
  }

  TRandom rndm;
  Long64_t nentries = t->GetEntriesFast();
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    t->GetEntry(jentry);
    if (fabs(vz) > 30) continue;
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
      
      for (int i = 0; i < smears; i++) {
        for (int j = 0; j < smears; j++) {
          for (int k = 0; k < smears; k++) {
            for (int l = 0; l < smears; l++) {
              pfrac = (pfrac_high-pfrac_low)/(float)smears*i + pfrac_low;
              efrac = (efrac_high-efrac_low)/(float)smears*j + efrac_low;
              econst = (econst_high-econst_low)/(float)smears*k + econst_low;
              escale = (escale_high-escale_low)/(float)smears*l + escale_low;
              
              float eta1fl = rndm.Gaus(eta1,pfrac/(e1)); 
              float phi1fl = rndm.Gaus(phi1,pfrac/(e1)); 
              float eta2fl = rndm.Gaus(eta2,pfrac/(e2)); 
              float phi2fl = rndm.Gaus(phi2,pfrac/(e2)); 

              float e1fl = rndm.Gaus(e1*escale,econst+efrac/TMath::Sqrt(e1));
              float e2fl = rndm.Gaus(e2*escale,econst+efrac/TMath::Sqrt(e2));
              float pt1fl = e1fl/cosh(eta1fl);
              float pt2fl = e2fl/cosh(eta2fl);

              if (eta1fl > etamax || eta1fl < etamin) continue;
              if (eta2fl > etamax || eta2fl < etamin) continue;

              TLorentzVector pho1fl;
              TLorentzVector pho2fl;
              pho1fl.SetPtEtaPhiE(pt1fl,eta1fl,phi1fl,e1fl);
              pho2fl.SetPtEtaPhiE(pt2fl,eta2fl,phi2fl,e2fl);
              TLorentzVector phofl = pho1fl + pho2fl;
              if(pt1fl < ptcut || pt2fl < ptcut) continue;
              if (phofl.Eta() > etamax || phofl.Eta() < etamin) continue;

              float energyimbal = abs(pho1fl.E() - pho2fl.E())/(pho1fl.E() + pho2fl.E());
              if(energyimbal > 0.6) continue;

              pt = phofl.Pt();

              int bin;// = (int) pt;
              if (pt < 2.5) bin = 0;
              else if (pt >= 2.5 && pt < 3.5) bin = 1;
              else if (pt >= 3.5) bin = 2;

              int abin;
              if (energyimbal < 0.08) abin = 0;
              else if (energyimbal >= 0.08 && energyimbal < 0.2) abin = 1;
              else if (energyimbal >= 0.2) abin = 2;

              float mass = phofl.M(); 

              hmass_pi0[bin][abin][i][j][k][l]->Fill(mass);
              hmass_eta[bin][abin][i][j][k][l]->Fill(mass);
            }
          }
        }
      }
    }
  }

  wf->cd();
  for(int ip=0; ip < nPtBins; ip++){
    for (int a =0; a < nAlphaBins; a++) {
      for(int i=0; i < smears; i++){
        for(int j=0; j < smears; j++){
          for(int k=0; k < smears; k++){
            for(int l=0; l < smears; l++){
              hmass_pi0[ip][a][i][j][k][l]->Write();
              hmass_eta[ip][a][i][j][k][l]->Write();
            }
          }
        }
      }
    }
  }
}
