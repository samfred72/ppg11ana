#ifndef CALOANA_H__
#define CALOANA_H__

// Utility
#include <vector>
#include <fstream>
#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH2.h>
#include <TGraph2D.h>
#include <TF2.h>
#include <cassert>
#include <sstream>
#include <string>
#include <TLorentzVector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos

// Fun4All
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

// Event
#include <Event/Event.h>
#include <Event/packet.h>
#include <ffaobjects/EventHeaderv1.h>

//Trigger
#include <calotrigger/TriggerRunInfov1.h>
#include <calotrigger/TriggerAnalyzer.h>
#include <calotrigger/LL1Out.h>
#include <calotrigger/LL1Outv1.h>
#include <calotrigger/TriggerPrimitive.h>
#include <calotrigger/TriggerPrimitivev1.h>
#include <calotrigger/TriggerPrimitiveContainer.h>
#include <calotrigger/TriggerPrimitiveContainerv1.h>
#include <calotrigger/TriggerDefs.h>

// Jet base
#include <jetbase/Jet.h>
#include <jetbase/Jetv1.h>
#include <jetbase/JetAlgo.h>
#include <jetbase/JetMap.h>
#include <jetbase/JetContainer.h>
#include <jetbackground/TowerBackgroundv1.h>

// Global vertex
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertexMapv1.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMapv1.h>
#include <globalvertex/GlobalVertex.h>

// G4
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>

// G4Cells includes
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>

// Tower includes
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfov1.h>
#include <calobase/TowerInfoContainerSimv1.h>
#include <calobase/TowerInfoSimv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfov2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfov3.h>
#include <calobase/TowerInfoContainerv4.h>
#include <calobase/TowerInfov4.h>
#include <calobase/TowerInfoDefs.h>

// MBD
#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtSimContainerV1.h>
#include <mbd/MbdPmtHit.h>
#include <mbd/MbdGeom.h>

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>

// phool
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#pragma GCC diagnostic pop

//#include <HepMC/GenVertex.h>             // for GenVertex, GenVertex::partic...
//#include <HepMC/GenParticle.h>           // for GenParticle
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

// Centrality MB
#include <centrality/CentralityInfo.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <centrality/CentralityInfov1.h>

#include <ffarawobjects/Gl1Packet.h>

// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;
class TH2F;
class TH1;
class Gl1Packet;
class PHG4Shower;
class TriggerPrimitive;
class TriggerPrimitiveContainer;
class LL1Out;


class CaloAna : public SubsysReco
{
 public:
  //! constructor
  CaloAna(const std::string &name = "CaloAna",  const char* outname="DST-00021615-0000.root", bool _isMC=false, bool _iHI=false);

  //! destructor
  virtual ~CaloAna();

  //! full initialization
  int Init(PHCompositeNode *);
  void InitOutputFile();
  void InitTree();
  
  //! event processing method
  int process_event(PHCompositeNode *);
  int ProcessGlobalEventInfo(PHCompositeNode *);
  void ProcessFillTruthParticle(PHCompositeNode *, float truthvz);

  //! end of run method
  int End(PHCompositeNode *);

  Double_t GetShiftedEta(float _vz, float _eta);

  int process_g4hits(PHCompositeNode *);
  int process_g4cells(PHCompositeNode *);
  int process_towers(PHCompositeNode *);
  int process_clusters(PHCompositeNode *);
  bool FindConversion(PHG4TruthInfoContainer *, int trackid, float energy);
  int Getpeaktime(TH1 *h);

  void SetMbdZVtxCut(float _zvtxcut)
  {
    if(_zvtxcut<=0) IsVtxCut = false;
    else IsVtxCut = true;
    vzcut = _zvtxcut;
  };

  void SetPdgId(int _pdgid)
  {
    PDGPID = _pdgid;
  };

  void SetMCFlags(bool _isPythia, bool _isSingleGun)
  {
    isPythia = _isPythia;
    isSingleGun = _isSingleGun;
    if(_isSingleGun) UseClusterTruthVertex = true;
  };
  
  void SetClusterVertexToZero(bool _isSetVtxZero)
  {
    UseClusterVertexZero = _isSetVtxZero;
  }

  void SetUseOfClusters(bool _doClusters)
  {
    doClusters = _doClusters;
  }

  void SetUseOfJets(bool _doJets)
  {
    doJets = _doJets;
  }

  void SetRunNumber(int RunNumber)
  {
    m_run_number = RunNumber;
  };

  void SetScaledowns(int scaledowns[])
  {
    for(int i=0; i<64;i++){
      m_scaledowns[i] = scaledowns[i];
      std::cout << "prescale set " << i << " : " << m_scaledowns[i] << std::endl;
    }
  };

  void SetTriggerEmulator(bool _useEmulator)
  {
    useEmulator = _useEmulator;
  };

  void Detector(const std::string &name) { detector = name; }
  void ProcessClearBranchVar();

 protected:
  std::string detector;
  std::string outfilename;
  Fun4AllHistoManager *hm = nullptr;
  TFile *outfile = nullptr;
  TTree *towerntuple = nullptr;
  TH1D * profilehist;
  TH1D * histlumicount;
  Long64_t mbdlivecount=0;

  int m_run_number=-999;
  int m_scaledowns[64];
  float m_emcal_energy[24576];
  int   m_emcal_etabin[24576];
  int   m_emcal_phibin[24576];
  float m_emcal_total_energy = 0;
  float m_mbd_energy_south = 0;
  float m_mbd_energy_north = 0;
  float _mbd_charge_threshold = 0.4;
  int m_mbd_hits_south = 0;
  int m_mbd_hits_north = 0;
  float m_mbd_charge_south[64];
  float m_mbd_charge_north[64];
  float m_mbd_time_south[64];
  float m_mbd_time_north[64];

  Short_t nClusters=0;
  Short_t nClusters_mother=0;
  int count=0;
  int count_truthvz=0;
  int count_mbdvtx=0;
  int count_mbdhit=0;
  Short_t nPairs=0;
  Short_t nPairsMother=0;
  static const int Max_diphoton_size = 10000;
  static const int Max_photon_size =   10000;
  Short_t idx_photon1[Max_diphoton_size];
  Short_t idx_photon2[Max_diphoton_size];
  Short_t idx_photon1_mother[Max_diphoton_size];
  Short_t idx_photon2_mother[Max_diphoton_size];
  short reco_matched_truth_idx[Max_diphoton_size];
  std::vector<float> photon_profvar;

  int m_UseAltZVertex = 1;
  float diphoton_energyimbal[Max_diphoton_size];
  float diphoton_energyimbal_mother[Max_diphoton_size];
  float diphoton_dR[Max_diphoton_size];
  float diphoton_dR_mother[Max_diphoton_size];
  float photon_prob[Max_photon_size];
  float photon_prob_merged_cluster[Max_photon_size];
  float photon_prob_mother[Max_photon_size];
  float photon_prob_mother_merged_cluster[Max_photon_size];
  std::vector<std::vector<float>> profvec;
  std::vector<std::vector<float>> nearbytowers;
  std::vector<std::vector<float>> nearbytowers_mother;
  float eta_cluster_center[Max_photon_size];
  float phi_cluster_center[Max_photon_size];
  float eta_cluster_mother_center[Max_photon_size];
  float phi_cluster_mother_center[Max_photon_size];
  float e1[Max_photon_size];
  float e2[Max_photon_size];
  float e3[Max_photon_size];
  float e4[Max_photon_size];
  float et1[Max_photon_size];
  float et2[Max_photon_size];
  float et3[Max_photon_size];
  float et4[Max_photon_size];
  float en[Max_photon_size];
  float rr[Max_photon_size];
  float ddy[Max_photon_size];
  float ddz[Max_photon_size];
  float theta[Max_photon_size];
  int iycg0[Max_photon_size];
  int izcg0[Max_photon_size];
  
  float mother_e1[Max_photon_size];
  float mother_e2[Max_photon_size];
  float mother_e3[Max_photon_size];
  float mother_e4[Max_photon_size];
  float mother_et1[Max_photon_size];
  float mother_et2[Max_photon_size];
  float mother_et3[Max_photon_size];
  float mother_et4[Max_photon_size];
  float mother_en[Max_photon_size];
  float mother_rr[Max_photon_size];
  float mother_ddy[Max_photon_size];
  float mother_ddz[Max_photon_size];
  float mother_theta[Max_photon_size];
  int mother_iycg0[Max_photon_size];
  int mother_izcg0[Max_photon_size];
  int nsubclusters[Max_photon_size];

  TClonesArray* truth_photon1_4mom;
  TClonesArray* truth_photon2_4mom;
  TClonesArray* truth_diphoton_4mom;
  TClonesArray* photon_4mom;
  TClonesArray* photon_4mom_ecore;
  TClonesArray* photon_4mom_mother;
  TClonesArray* photon_4mom_mother_ecore;
  TClonesArray* diphoton_4mom;
  TClonesArray* diphoton_4mom_mother;
  TClonesArray* photon1_4mom_fit;
  TClonesArray* photon2_4mom_fit;
  TClonesArray* diphoton_4mom_fit;
  float fit_chisq[Max_photon_size];
  float fitres[Max_photon_size];

  float vx, vy, vz;
  float vzcut;
  bool IsVtxCut = true;
  bool ScaledTriggerBit[64];
  bool LiveTriggerBit[64];
  long long int count_raw[64];
  long long int count_live[64];
  long long int count_scaled[64]; 

  float _emcal_cluster_emincut = 0.5;
  float diphotonmasslow = 0.;
  float diphotonmasshigh = 2.5;

  bool isMC;
  bool isHI;
  int npart;
  int ncoll;
  float bimp;
  int centbin;
  bool isMinBias;

  int Centrality;
  float Cent_impactparam;

  int processId;
  static const int nParticleTruth = 10000;
  float truth_vz=-999;
  float truth_vx=-999;
  float truth_vy=-999;
  int truthpar_n = 0;
  int truth_id[nParticleTruth];
  bool isconverted[nParticleTruth];
  bool isconverted1[nParticleTruth];
  bool isconverted2[nParticleTruth];
  int truth_track_id[nParticleTruth];
  bool found_decay1[nParticleTruth];
  bool found_decay2[nParticleTruth];

  int PDGPID = 111;
  TLorentzVector truthphoton1;
  TLorentzVector truthphoton2;

  bool doClusters =false;
  bool doJets =false;
  bool doEMCal = false;
  bool doHCal = false;
  bool doMBD = false;
  bool doZDC = false;
  bool isPythia = false;
  bool UseClusterTruthVertex = false;
  bool UseClusterVertexZero = false;
  bool isSingleGun = false;

  TriggerAnalyzer *triggeranalyzer{nullptr};
  bool useEmulator = false;

  const int inputDimx{7};
  const int inputDimy{7};
    
  float _jet_minetcut = 5;
  float _jet_etacut= 2;

  static const int Max_jet_size = 1000;
  int nJets=0;
  float jet_e[Max_jet_size];
  float jet_pt[Max_jet_size];
  float jet_eta[Max_jet_size];
  float jet_phi[Max_jet_size];
  float m_jetfrac_ihcal[Max_jet_size];
  float m_jetfrac_ohcal[Max_jet_size];
  float m_jetfrac_emcal[Max_jet_size];

  int nTruthJets=0;
  float truth_jet_e[Max_jet_size];
  float truth_jet_pt[Max_jet_size];
  float truth_jet_eta[Max_jet_size];
  float truth_jet_phi[Max_jet_size];

};

#endif
