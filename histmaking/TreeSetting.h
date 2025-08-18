#ifndef TREESETTING_h
#define TREESETTING_h

#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include <TLorentzVector.h>
#include <vector>

static const long int MAXCLUSTERS=1000;
Int_t           RunNumber;
Float_t         vz;
Float_t         emcal_total_energy;
TClonesArray    *truth_photon1_4mom;
TClonesArray    *truth_photon2_4mom;
TClonesArray    *truth_diphoton_4mom;
Float_t         truth_vz;
Int_t           truthpar_n;
Int_t           truth_id[MAXCLUSTERS];   //[truthpar_n]
Bool_t          truth_isconverted[MAXCLUSTERS];   //[truthpar_n]
Bool_t          truth_isconverted1[MAXCLUSTERS];   //[truthpar_n]
Bool_t          truth_isconverted2[MAXCLUSTERS];   //[truthpar_n]
Bool_t          truth_found_decay1[MAXCLUSTERS];   //[truthpar_n]
Bool_t          truth_found_decay2[MAXCLUSTERS];   //[truthpar_n]
Short_t         nClusters;
Short_t         nClusters_mother;
TClonesArray    *photon_4mom;
TClonesArray    *photon_4mom_mother;
Float_t         photon_prob[MAXCLUSTERS];   //[nClusters]
Float_t         photon_prob_merged_cluster[MAXCLUSTERS];   //[nClusters]
Float_t         photon_prob_mother[MAXCLUSTERS];   //[nClusters_mother]
Float_t         photon_prob_mother_merged_cluster[MAXCLUSTERS];   //[nClusters_mother]
Float_t         eta_cluster_mother_center[MAXCLUSTERS];   //[nClusters_mother]
Float_t         phi_cluster_mother_center[MAXCLUSTERS];   //[nClusters_mother]
Float_t         eta_cluster_center[MAXCLUSTERS];   //[nClusters_mother]
Float_t         phi_cluster_center[MAXCLUSTERS];   //[nClusters_mother]
std::vector<std::vector<float> > *nearbytowers;
std::vector<std::vector<float> > *nearbytowers_mother;
Short_t         nPairs;
TClonesArray    *diphoton_4mom;
Float_t         diphoton_energyimbal[MAXCLUSTERS];   //[nPairs]
Float_t         diphoton_dR[MAXCLUSTERS];   //[nPairs]
Short_t         idx_photon1[MAXCLUSTERS];   //[nPairs]
Short_t         idx_photon2[MAXCLUSTERS];   //[nPairs]
Short_t         nPairsMother;
TClonesArray    *diphoton_4mom_mother;
Short_t         idx_photon1_mother[MAXCLUSTERS];   //[nPairsMother]
Short_t         idx_photon2_mother[MAXCLUSTERS];   //[nPairsMother]
Short_t         reco_matched_truth_idx[MAXCLUSTERS];   //[nPairs]

TClonesArray    *diphoton_4mom_fit;
TClonesArray    *photon1_4mom_fit;
TClonesArray    *photon2_4mom_fit;
Float_t         fit_chisq[MAXCLUSTERS];
Bool_t          ScaledTriggerBit[64];
Bool_t          LiveTriggerBit[64];

// List of branches
TBranch        *b_RunNumber;   //!
TBranch        *b_vz;   //!
TBranch        *b_emcal_total_energy;   //!
TBranch        *b_truth_photon1_4mom;   //!
TBranch        *b_truth_photon2_4mom;   //!
TBranch        *b_truth_diphoton_4mom;   //!
TBranch        *b_truth_vz;   //!
TBranch        *b_truthpar_n;   //!
TBranch        *b_truth_id;   //!
TBranch        *b_truth_isconverted;   //!
TBranch        *b_truth_isconverted1;   //!
TBranch        *b_truth_isconverted2;   //!
TBranch        *b_truth_found_decay1;   //!
TBranch        *b_truth_found_decay2;   //!
TBranch        *b_nClusters;   //!
TBranch        *b_nClusters_mother;   //!
TBranch        *b_photon_4mom;   //!
TBranch        *b_photon_4mom_mother;   //!
TBranch        *b_photon_prob_mother;   //!
TBranch        *b_photon_prob_mother_merged_cluster;   //!
TBranch        *b_photon_prob;   //!
TBranch        *b_photon_prob_merged_cluster;   //!
TBranch        *b_eta_cluster_center;   //!
TBranch        *b_phi_cluster_center;   //!
TBranch        *b_eta_cluster_mother_center;   //!
TBranch        *b_phi_cluster_mother_center;   //!
TBranch        *b_nearbytowers;   //!
TBranch        *b_nearbytowers_mother;   //!
TBranch        *b_nPairs;   //!
TBranch        *b_diphoton_4mom;   //!
TBranch        *b_diphoton_energyimbal;   //!
TBranch        *b_diphoton_dR;   //!
TBranch        *b_idx_photon1;   //!
TBranch        *b_idx_photon2;   //!
TBranch        *b_nPairsMother;   //!
TBranch        *b_diphoton_4mom_mother;   //!
TBranch        *b_idx_photon1_mother;   //!
TBranch        *b_idx_photon2_mother;   //!
TBranch        *b_reco_matched_truth_idx;   //!

TBranch        *b_diphoton_4mom_fit;
TBranch        *b_photon1_4mom_fit;
TBranch        *b_photon2_4mom_fit;
TBranch        *b_fit_chisq;
TBranch        *b_ScaledTriggerBit;
TBranch        *b_LiveTriggerBit;

void treesetup(TTree* tree, bool isMC = 0){ 
  
  if (!tree){
    std::cerr << "ERROR!! TTree does not exist...." << std::endl;
    return;
  }

  truth_photon1_4mom = 0;
  truth_photon2_4mom = 0;
  truth_diphoton_4mom = 0;
  photon_4mom = 0;
  photon_4mom_mother = 0;
  nearbytowers = 0;
  nearbytowers_mother = 0;
  diphoton_4mom = 0;
  diphoton_4mom_mother = 0;

  tree->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
  tree->SetBranchAddress("vz", &vz, &b_vz);
  tree->SetBranchAddress("emcal_total_energy", &emcal_total_energy, &b_emcal_total_energy);
  if(isMC){
    tree->SetBranchAddress("truth_photon1_4mom", &truth_photon1_4mom, &b_truth_photon1_4mom);
    tree->SetBranchAddress("truth_photon2_4mom", &truth_photon2_4mom, &b_truth_photon2_4mom);
    tree->SetBranchAddress("truth_diphoton_4mom", &truth_diphoton_4mom, &b_truth_diphoton_4mom);
    tree->SetBranchAddress("truth_vz", &truth_vz, &b_truth_vz);
    tree->SetBranchAddress("truthpar_n", &truthpar_n, &b_truthpar_n);
    tree->SetBranchAddress("truth_id", truth_id, &b_truth_id);
    tree->SetBranchAddress("truth_isconverted", truth_isconverted, &b_truth_isconverted);
    tree->SetBranchAddress("truth_isconverted1", truth_isconverted1, &b_truth_isconverted1);
    tree->SetBranchAddress("truth_isconverted2", truth_isconverted2, &b_truth_isconverted2);
    tree->SetBranchAddress("truth_found_decay1", truth_found_decay1, &b_truth_found_decay1);
    tree->SetBranchAddress("truth_found_decay2", truth_found_decay2, &b_truth_found_decay2);
    tree->SetBranchAddress("reco_matched_truth_idx", reco_matched_truth_idx, &b_reco_matched_truth_idx);
  }
  tree->SetBranchAddress("nClusters", &nClusters, &b_nClusters);
  tree->SetBranchAddress("nClusters_mother", &nClusters_mother, &b_nClusters_mother);
  tree->SetBranchAddress("photon_4mom", &photon_4mom, &b_photon_4mom);
  tree->SetBranchAddress("photon_4mom_mother", &photon_4mom_mother, &b_photon_4mom_mother);
  tree->SetBranchAddress("photon_prob_mother", &photon_prob_mother, &b_photon_prob_mother);
  tree->SetBranchAddress("photon_prob_mother_merged_cluster", &photon_prob_mother_merged_cluster, &b_photon_prob_mother_merged_cluster);
  tree->SetBranchAddress("photon_prob", photon_prob, &b_photon_prob);
  tree->SetBranchAddress("photon_prob_merged_cluster", photon_prob_merged_cluster, &b_photon_prob_merged_cluster);
  tree->SetBranchAddress("eta_cluster_center", eta_cluster_center, &b_eta_cluster_center);
  tree->SetBranchAddress("phi_cluster_center", phi_cluster_center, &b_phi_cluster_center);
  tree->SetBranchAddress("eta_cluster_mother_center", &eta_cluster_mother_center, &b_eta_cluster_mother_center);
  tree->SetBranchAddress("phi_cluster_mother_center", &phi_cluster_mother_center, &b_phi_cluster_mother_center);
  tree->SetBranchAddress("nearbytowers", &nearbytowers, &b_nearbytowers);
  tree->SetBranchAddress("nearbytowers_mother", &nearbytowers_mother, &b_nearbytowers_mother);
  tree->SetBranchAddress("nPairs", &nPairs, &b_nPairs);
  tree->SetBranchAddress("diphoton_4mom", &diphoton_4mom, &b_diphoton_4mom);
  tree->SetBranchAddress("diphoton_energyimbal", diphoton_energyimbal, &b_diphoton_energyimbal);
  tree->SetBranchAddress("diphoton_dR", diphoton_dR, &b_diphoton_dR);
  tree->SetBranchAddress("idx_photon1", idx_photon1, &b_idx_photon1);
  tree->SetBranchAddress("idx_photon2", idx_photon2, &b_idx_photon2);
  tree->SetBranchAddress("nPairsMother", &nPairsMother, &b_nPairsMother);
  tree->SetBranchAddress("diphoton_4mom_mother", &diphoton_4mom_mother, &b_diphoton_4mom_mother);
  tree->SetBranchAddress("idx_photon1_mother", idx_photon1_mother, &b_idx_photon1_mother);
  tree->SetBranchAddress("idx_photon2_mother", idx_photon2_mother, &b_idx_photon2_mother);
  tree->SetBranchAddress("diphoton_4mom_fit", &diphoton_4mom_fit, &b_diphoton_4mom_fit);
  tree->SetBranchAddress("photon1_4mom_fit", &photon1_4mom_fit, &b_photon1_4mom_fit);
  tree->SetBranchAddress("photon2_4mom_fit", &photon2_4mom_fit, &b_photon2_4mom_fit);
  tree->SetBranchAddress("fit_chisq", fit_chisq, &b_fit_chisq);
  tree->SetBranchAddress("ScaledTriggerBit", ScaledTriggerBit, &b_ScaledTriggerBit);
  tree->SetBranchAddress("LiveTriggerBit", LiveTriggerBit, &b_LiveTriggerBit);
}

#endif 
