#ifndef EVENTFILTER_H__
#define EVENTFILTER_H__

// Utility
#include <vector>
#include <fstream>
#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH2.h>
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

// Jet base
#include <jetbase/Jet.h>
#include <jetbase/Jetv1.h>
#include <jetbase/JetAlgo.h>
#include <jetbase/JetMap.h>
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
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>

// phool
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

// HepMC
#include <phhepmc/PHHepMCGenEventMap.h>


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#pragma GCC diagnostic pop

//#include <HepMC/GenVertex.h>             // for GenVertex, GenVertex::partic...
//#include <HepMC/GenParticle.h>           // for GenParticle
#include <phhepmc/PHHepMCGenEvent.h>
//#include <phhepmc/PHHepMCGenEventMap.h>

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

class EventFilter : public SubsysReco
{
 public:
  //! constructor
  EventFilter(const std::string &name = "EventFilter",  const char* outname="DST-00021615-0000.root", bool _isMC=false, bool _iHI=false);

  //! destructor
  virtual ~EventFilter();

  //! full initialization
  int Init(PHCompositeNode *);
  void InitOutputFile();
  void InitTree();
  
  //! event processing method
  int process_event(PHCompositeNode *);
  int ProcessGlobalEventInfo(PHCompositeNode *);
  void ProcessFillTruthParticle(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  int process_g4hits(PHCompositeNode *);
  int process_g4cells(PHCompositeNode *);
  int process_towers(PHCompositeNode *);
  int process_clusters(PHCompositeNode *);
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

  int m_run_number=-999;
  int m_scaledowns[64];
  float m_emcal_energy[24576];
  float m_emcal_total_energy = 0;
  float m_mbd_energy_south = 0;
  float m_mbd_energy_north = 0;
  float _mbd_charge_threshold = 0.5;
  int m_mbd_hits_south = 0;
  int m_mbd_hits_north = 0;
  float m_mbd_charge_south[64];
  float m_mbd_charge_north[64];
  float m_mbd_time_south[64];
  float m_mbd_time_north[64];

  Short_t nClusters=0;
  int count=0;
  int count_mbdvtx=0;
  Short_t nPairs=0;
  static const int Max_diphoton_size = 10000;
  static const int Max_photon_size = 10000;
  Short_t idx_photon1[Max_diphoton_size];
  Short_t idx_photon2[Max_diphoton_size];
  float reco_matched_truth_idx[Max_diphoton_size];

  int m_UseAltZVertex = 1;
  float diphoton_energyimbal[Max_diphoton_size];
  float diphoton_dR[Max_diphoton_size];
  float photon_prob[Max_photon_size];

  TClonesArray* photon_4mom;
  TClonesArray* diphoton_4mom;

  float vx, vy, vz;
  float vzcut;
  bool IsVtxCut = true;
  bool ScaledTriggerBit[64];
  bool LiveTriggerBit[64];
  long long int count_raw[64];
  long long int count_live[64];
  long long int count_scaled[64]; 

  float _emcal_cluster_emincut = 0.5;
  float diphotonmasslow = 0.01;
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
  float truth_energy[nParticleTruth];
  float truth_pt[nParticleTruth];
  float truth_eta[nParticleTruth];
  float truth_phi[nParticleTruth];
  int truth_id[nParticleTruth];

  int PDGPID = 22;
  TLorentzVector truthphoton1;
  TLorentzVector truthphoton2;
  bool found_decay;

  bool doClusters = true;
  bool doEMCal = false;
  bool doHCal = false;
  bool doMBD = false;
  bool doZDC = false;
  bool isPythia = false;
  bool UseClusterTruthVertex = false;
  bool isSingleGun = false;

  TriggerAnalyzer *triggeranalyzer{nullptr};
  bool useEmulator = false;
};

#endif
