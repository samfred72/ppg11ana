#include <iostream>
#include "CaloAna.h"
using namespace std;

  CaloAna::CaloAna(const std::string& name, const char* outname, bool _isMC, bool _isHI)
: SubsysReco(name)
  , detector("HCALIN")
{
  outfilename = Form("%s",outname);
  isMC = _isMC;
  isHI = _isHI;
}

CaloAna::~CaloAna()
{
    //delete hm;
    delete towerntuple;
}

int CaloAna::Init(PHCompositeNode*)
{ 
  triggeranalyzer = new TriggerAnalyzer();
  triggeranalyzer->UseEmulator(useEmulator);
  try {
    InitOutputFile();
    InitTree();
  }  
  catch (const std::exception& e){
    std::cerr << "CaloAna::Init - Exception during init! For " << e.what() << " return -1" << std::endl;
    return -1;
  }
  TFile * fhist = TFile::Open("/sphenix/user/samfred/ppg11ana/Profile/newprof/profilehists.root","READ");
  profilehist = (TH1D*)fhist->Get("energyprof1D");

  return 0;
}

void CaloAna::InitOutputFile(){
  std::cout << "output filename : " << outfilename.c_str() << std::endl;
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  if(!outfile || outfile->IsZombie()){
    throw std::runtime_error("Failed open file");
  }
}

void CaloAna::InitTree(){
    towerntuple = new TTree("towerntup", "Ntuple");
    towerntuple->Branch("RunNumber",&m_run_number);
    towerntuple->Branch("vz",&vz);
    if(doEMCal) towerntuple->Branch("emcal_total_energy",&m_emcal_total_energy);
    if(!isSingleGun){
      towerntuple->Branch("mbd_total_charge_south",&m_mbd_energy_south);
      towerntuple->Branch("mbd_total_charge_north",&m_mbd_energy_north);
      towerntuple->Branch("mbd_nhits_south",&m_mbd_hits_south);
      towerntuple->Branch("mbd_nhits_north",&m_mbd_hits_north);
    }
    if(doMBD){
      towerntuple->Branch("mbd_charge_south",m_mbd_charge_south);
      towerntuple->Branch("mbd_charge_north",m_mbd_charge_north);
      towerntuple->Branch("mbd_time_south",m_mbd_time_south);
      towerntuple->Branch("mbd_time_north",m_mbd_time_north);
    }

    if(isMC){
      truth_photon1_4mom= new TClonesArray("TLorentzVector", Max_photon_size);
      truth_photon2_4mom= new TClonesArray("TLorentzVector", Max_photon_size);
      truth_diphoton_4mom= new TClonesArray("TLorentzVector", Max_photon_size);
      towerntuple->Branch("truth_photon1_4mom", "TClonesArray", &truth_photon1_4mom, 32000, 0); 
      towerntuple->Branch("truth_photon2_4mom", "TClonesArray", &truth_photon2_4mom, 32000, 0); 
      towerntuple->Branch("truth_diphoton_4mom", "TClonesArray", &truth_diphoton_4mom, 32000, 0); 
      towerntuple->Branch("truth_vz",&truth_vz);
      towerntuple->Branch("truthpar_n",&truthpar_n);
      towerntuple->Branch("truth_id",truth_id,"truth_id[truthpar_n]/I");
      towerntuple->Branch("truth_isconverted",isconverted,"truth_isconverted[truthpar_n]/O");
      towerntuple->Branch("truth_found_decay1",found_decay1,"truth_found_decay1[truthpar_n]/O");
      towerntuple->Branch("truth_found_decay2",found_decay2,"truth_found_decay2[truthpar_n]/O");
      towerntuple->Branch("truth_isconverted1",isconverted1,"truth_isconverted1[truthpar_n]/O");
      towerntuple->Branch("truth_isconverted2",isconverted2,"truth_isconverted2[truthpar_n]/O");
      if(isPythia) towerntuple->Branch("processId",&processId);
    }
    else if(!isMC){
      towerntuple->Branch("ScaledTriggerBit",ScaledTriggerBit,"ScaledTriggerBit[64]/O");
      towerntuple->Branch("LiveTriggerBit",LiveTriggerBit,"LiveTriggerBit[64]/O");
      towerntuple->Branch("Scaledowns",m_scaledowns,"Scaledowns[64]/I");
      towerntuple->Branch("count_raw",count_raw,"count_raw[64]/L");
      towerntuple->Branch("count_live",count_live,"count_live[64]/L");
      towerntuple->Branch("count_scaled",count_scaled,"count_scaled[64]/L");
      histlumicount = new TH1D("histlumicount",";;",1,0,1);
    }


    if(doClusters){
      photon_4mom= new TClonesArray("TLorentzVector", Max_photon_size);
      photon_4mom_mother= new TClonesArray("TLorentzVector", Max_photon_size);
      towerntuple->Branch("nClusters",&nClusters, "nClusters/S");
      towerntuple->Branch("nClusters_mother",&nClusters_mother, "nClusters_mother/S");
      towerntuple->Branch("photon_4mom", "TClonesArray", &photon_4mom, 32000, 0); 
      towerntuple->Branch("photon_4mom_mother", "TClonesArray", &photon_4mom_mother, 32000, 0); // used to be 32000
      towerntuple->Branch("photon_prob_mother",photon_prob_mother,"photon_prob_mother[nClusters_mother]/F");
      towerntuple->Branch("photon_prob_mother_merged_cluster",photon_prob_mother_merged_cluster,"photon_prob_mother_merged_cluster[nClusters_mother]/F");
      towerntuple->Branch("photon_prob",photon_prob,"photon_prob[nClusters]/F");
      towerntuple->Branch("photon_prob_merged_cluster",photon_prob_merged_cluster,"photon_prob_merged_cluster[nClusters]/F");
      towerntuple->Branch("eta_cluster_center",eta_cluster_center,"eta_cluster_center[nClusters]/F");
      towerntuple->Branch("phi_cluster_center",phi_cluster_center,"phi_cluster_center[nClusters]/F");
      towerntuple->Branch("eta_cluster_mother_center",eta_cluster_mother_center,"eta_cluster_mother_center[nClusters_mother]/F");
      towerntuple->Branch("phi_cluster_mother_center",phi_cluster_mother_center,"phi_cluster_mother_center[nClusters_mother]/F");
      towerntuple->Branch("nearbytowers",&nearbytowers);
      towerntuple->Branch("nearbytowers_mother",&nearbytowers_mother);
      diphoton_4mom = new TClonesArray("TLorentzVector", Max_diphoton_size);
      towerntuple->Branch("nPairs",&nPairs,"nPairs/S");
      towerntuple->Branch("diphoton_4mom", "TClonesArray", &diphoton_4mom, 32000, 0); 
      towerntuple->Branch("diphoton_energyimbal",diphoton_energyimbal,"diphoton_energyimbal[nPairs]/F");
      towerntuple->Branch("diphoton_dR",diphoton_dR,"diphoton_dR[nPairs]/F");
      towerntuple->Branch("idx_photon1",idx_photon1,"idx_photon1[nPairs]/S");
      towerntuple->Branch("idx_photon2",idx_photon2,"idx_photon2[nPairs]/S");
      
      diphoton_4mom_mother = new TClonesArray("TLorentzVector", Max_diphoton_size);
      diphoton_4mom_fit = new TClonesArray("TLorentzVector", Max_diphoton_size);
      photon1_4mom_fit = new TClonesArray("TLorentzVector", Max_diphoton_size);
      photon2_4mom_fit = new TClonesArray("TLorentzVector", Max_diphoton_size);
      towerntuple->Branch("nPairsMother",&nPairsMother,"nPairsMother/S");
      towerntuple->Branch("diphoton_4mom_mother", "TClonesArray", &diphoton_4mom_mother, 32000, 0); // used to be 32000
      towerntuple->Branch("diphoton_energyimbal_mother",diphoton_energyimbal_mother,"diphoton_energyimbal_mother[nPairsMother]/F");
      towerntuple->Branch("diphoton_dR_mother",diphoton_dR_mother,"diphoton_dR_mother[nPairsMother]/F");
      towerntuple->Branch("idx_photon1_mother",idx_photon1_mother,"idx_photon1_mother[nPairsMother]/S");
      towerntuple->Branch("idx_photon2_mother",idx_photon2_mother,"idx_photon2_mother[nPairsMother]/S");
      towerntuple->Branch("photon1_4mom_fit","TClonesArray", &photon1_4mom_fit, 32000, 0);
      towerntuple->Branch("photon2_4mom_fit","TClonesArray", &photon2_4mom_fit, 32000, 0);
      towerntuple->Branch("diphoton_4mom_fit","TClonesArray", &diphoton_4mom_fit, 32000, 0);
      towerntuple->Branch("fit_chisq",&fit_chisq,"fit_chisq[nPairsMother]/F");
      towerntuple->Branch("fitres",&fitres,"fitres[nPairsMother]/F");
      if(isMC){
        towerntuple->Branch("reco_matched_truth_idx",reco_matched_truth_idx,"reco_matched_truth_idx[nPairs]/S");
      }
    }
    
    if(doJets){
      towerntuple->Branch("nJets",&nJets,"nJets/I");
      towerntuple->Branch("jet_e",jet_e,"jet_e[nJets]/F");
      towerntuple->Branch("jet_pt",jet_pt,"jet_pt[nJets]/F");
      towerntuple->Branch("jet_eta",jet_eta,"jet_eta[nJets]/F");
      towerntuple->Branch("jet_phi",jet_phi,"jet_phi[nJets]/F");
      towerntuple->Branch("jetfrac_ihcal",m_jetfrac_ihcal,"jetfrac_ihcal[nJets]/F");
      towerntuple->Branch("jetfrac_ohcal",m_jetfrac_ohcal,"jetfrac_ohcal[nJets]/F");
      towerntuple->Branch("jetfrac_emcal",m_jetfrac_emcal,"jetfrac_emcal[nJets]/F");
      if(isMC){
        towerntuple->Branch("nTruthJets",&nTruthJets,"nTruthJets/I");
        towerntuple->Branch("truth_jet_e",truth_jet_e,"truth_jet_e[nTruthJets]/F");
        towerntuple->Branch("truth_jet_pt",truth_jet_pt,"truth_jet_pt[nTruthJets]/F");
        towerntuple->Branch("truth_jet_eta",truth_jet_eta,"truth_jet_eta[nTruthJets]/F");
        towerntuple->Branch("truth_jet_phi",truth_jet_phi,"truth_jet_phi[nTruthJets]/F");
      }
    }

}

int CaloAna::process_event(PHCompositeNode* topNode)
{
    if(!topNode){
      std::cerr << "CaloAna::Init - topnode PHCompositeNode not valid! return -1" << std::endl;
      return -1;
    }
    if(count % 100 == 0) std::cout << "event : " << count  << std::endl;
    count++;
    process_towers(topNode);
    return Fun4AllReturnCodes::EVENT_OK;
}

int CaloAna::process_towers(PHCompositeNode* topNode)
{
  if(!topNode){
    std::cerr << "CaloAna::process_towers - topnode not valid! return -1" << std::endl;
    return -1;
  }
  
  ProcessGlobalEventInfo(topNode);
  if(isMC){
    if(isSingleGun){
    //  if(fabs(truth_vz) > 50) return Fun4AllReturnCodes::ABORTEVENT;
    }
    //if(fabs(truth_vz) > 150) return Fun4AllReturnCodes::ABORTEVENT;
    count_truthvz++;
    ProcessFillTruthParticle(topNode,truth_vz);
//    if(truthpar_n ==0) return Fun4AllReturnCodes::ABORTEVENT;
  }
  
  if(IsVtxCut){
    if(fabs(vz)>vzcut) return Fun4AllReturnCodes::ABORTEVENT;
  }

  
  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  if (!geomEM) throw std::runtime_error("failed to find EMCAL TOWERGEOM node");

  // EMCAL energy
  TowerInfoContainer* offlinetowers = (isMC) ? static_cast<TowerInfoContainer*>(findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC")) : static_cast<TowerInfoContainer*>(findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC"));
  if(!offlinetowers) throw std::runtime_error("failed to find EMCAL TowerInfoContainer node");

  if(doEMCal){
    int size = offlinetowers->size();
    for (int channel = 0; channel < size;channel++)
    {
      TowerInfo* offlinetower = offlinetowers->get_tower_at_channel(channel);
      m_emcal_energy[channel] = offlinetower->get_energy();
      m_emcal_total_energy += m_emcal_energy[channel]; 
      if(offlinetower->get_isSaturated() ){
      std::cout << "energy of tower " << m_emcal_energy[channel]  << " and saturated flag : "  << offlinetower->get_isSaturated() << std::endl;
      }
      unsigned int towerkey = offlinetowers->encode_key(channel);
      int ieta = offlinetowers->getTowerEtaBin(towerkey);
      int iphi = offlinetowers->getTowerPhiBin(towerkey);

      m_emcal_etabin[channel] = ieta;
      m_emcal_phibin[channel] = iphi;
    }
  }

  // MBD energy and hits
  if(!isSingleGun){
    m_mbd_hits_south = 0;
    m_mbd_hits_north = 0;
    MbdPmtContainer * mbdpmts = (isMC) ? static_cast<MbdPmtContainer*>(findNode::getClass<MbdPmtSimContainerV1>(topNode, "MbdPmtContainer")) : static_cast<MbdPmtContainer*>(findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer"));
    if(mbdpmts){
      for (int ipmt=0; ipmt < 128; ipmt++){
        Float_t q = mbdpmts->get_pmt(ipmt)->get_q();
        Float_t t = mbdpmts->get_pmt(ipmt)->get_time();
        if(ipmt < 64) {
          m_mbd_energy_south += q;
          m_mbd_charge_south[ipmt] = q;
          m_mbd_time_south[ipmt] = t;
          if(q > _mbd_charge_threshold) {
            m_mbd_hits_south++;
          }
        }
        else if(ipmt >= 64){
          m_mbd_energy_north += q;
          m_mbd_charge_north[ipmt-64] = q;
          m_mbd_time_north[ipmt-64] = t;
          if(q > _mbd_charge_threshold) {
            m_mbd_hits_north++;
          }
        }
      }
    }
    if(m_mbd_hits_north>0 && m_mbd_hits_south>0)count_mbdhit++; 
  }

  if(doClusters)
  {
    RawClusterContainer *clustersEM = static_cast<RawClusterContainer*>(findNode::getClass<RawClusterContainer>(topNode, "CEMC_CLUSTERINFO"));
    if(!clustersEM) return Fun4AllReturnCodes::ABORTEVENT;
    RawClusterContainer *clustersEM_mother = static_cast<RawClusterContainer*>(findNode::getClass<RawClusterContainer>(topNode, "CEMC_CLUSTERINFO_MOTHER"));
    if(!clustersEM_mother) return Fun4AllReturnCodes::ABORTEVENT;
    RawClusterContainer::ConstIterator hiter;
    RawClusterContainer::ConstIterator hiter_piprob;
    RawClusterContainer::ConstIterator hiter2;
    RawClusterContainer::ConstIterator hiter2_piprob;
    RawClusterContainer::ConstIterator hiter_mother;
    RawClusterContainer::ConstIterator hiter_mother_piprob;
    RawClusterContainer::ConstIterator hiter2_mother;
    RawClusterContainer::ConstIterator hiter2_mother_piprob;
    RawClusterContainer::ConstRange begin_end = clustersEM->getClusters(); 
    RawClusterContainer::ConstRange begin_end_mother = clustersEM_mother->getClusters(); 

    nPairsMother=0;
    int _hiter1_mother=0;
    short photon1_idx_mother = 0;
    nClusters_mother=0;
    
    
    for(hiter_mother = begin_end_mother.first; hiter_mother != begin_end_mother.second; ++hiter_mother)
    {
      if (vz == -999 && isMC && isPythia) continue;
      _hiter1_mother++;
      RawCluster* cluster = hiter_mother->second;
      CLHEP::Hep3Vector vertex;
      if(isMC && UseClusterTruthVertex) vertex.set(truth_vx,truth_vy,truth_vz);
      else if(UseClusterVertexZero) vertex.set(0,0,0);
      else vertex.set(vx,vy,vz);
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetEVec(*cluster, vertex);
      
      float clPt = E_vec_cluster.perp();
      float clEta = E_vec_cluster.pseudoRapidity();
      float clPhi = E_vec_cluster.phi();
      float clE = E_vec_cluster.mag();
      //if(clPt < _emcal_cluster_emincut) continue;
      
      TLorentzVector vphoton1_mother;
      vphoton1_mother.SetPtEtaPhiE(clPt,clEta,clPhi,clE);
      
      photon1_idx_mother++;
      short photon2_idx_mother = 0;
      int _hiter2_mother=0;
      for(hiter2_mother = begin_end_mother.first; hiter2_mother != begin_end_mother.second; ++hiter2_mother){
        _hiter2_mother++;
        if(hiter_mother == hiter2_mother) continue;
        if(_hiter1_mother >= _hiter2_mother) continue;

        RawCluster* cluster2 = hiter2_mother->second;
        CLHEP::Hep3Vector E_vec_cluster2 = RawClusterUtility::GetEVec(*cluster2, vertex);
        float clPt2 = E_vec_cluster2.perp();
        float clEta2 = E_vec_cluster2.pseudoRapidity();
        float clPhi2 = E_vec_cluster2.phi();
        float clE2 = E_vec_cluster2.mag();

        if( abs(clEta2-clEta)<1e-4 && abs(clPhi2-clPhi)<1e-4) continue;
        //if (clPt2 < _emcal_cluster_emincut) continue;

        TLorentzVector vphoton2_mother;
        vphoton2_mother.SetPtEtaPhiE(clPt2,clEta2,clPhi2,clE2);

        TLorentzVector vdiphotonMother = vphoton1_mother+ vphoton2_mother;
        new ((*diphoton_4mom_mother)[nPairsMother]) TLorentzVector(vdiphotonMother);

        if(vdiphotonMother.M() < diphotonmasslow || vdiphotonMother.M()>diphotonmasshigh) continue;

        diphoton_energyimbal_mother[nPairsMother] = fabs(clE - clE2) / (clE + clE2);
        diphoton_dR_mother[nPairsMother] = vphoton1_mother.DeltaR(vphoton2_mother);
        
        photon2_idx_mother++;
        idx_photon1_mother[nPairsMother] = photon1_idx_mother -1;
        idx_photon2_mother[nPairsMother] = photon1_idx_mother -1 + photon2_idx_mother;
        nPairsMother++;
      }
        
      new ((*photon_4mom_mother)[nClusters_mother]) TLorentzVector(vphoton1_mother);
      photon_prob_mother[nClusters_mother] = cluster->get_prob();
      photon_prob_mother_merged_cluster[nClusters_mother] = cluster->get_merged_cluster_prob();

      std::vector<float> showershape = cluster->get_shower_shapes(0.07);
      eta_cluster_mother_center[nClusters_mother] = showershape.at(4);
      phi_cluster_mother_center[nClusters_mother] = showershape.at(5);

      int etacenter = std::floor(eta_cluster_mother_center[nClusters_mother] + 0.5);
      int phicenter = std::floor(phi_cluster_mother_center[nClusters_mother] + 0.5);
      
      std::vector<float> input;
      int vectorSize = inputDimx * inputDimy;
      input.resize(vectorSize, 0);

      int xlength = int((inputDimx - 1) / 2);
      int ylength = int((inputDimy - 1) / 2);

      // Graph for fitting
      float emax1 = 0;
      float emax2 = 0;
      TGraph2D * g_forfit = new TGraph2D();
      TH2D * h_cluster = new TH2D("h_cluster",";ieta;iphi",7,-0.5,6.5,7,-0.5,6.5);
      h_cluster->Sumw2();
      for (int ieta = etacenter - ylength; ieta <= etacenter + ylength; ieta++)
      {
        for (int iphi = phicenter - xlength; iphi <= phicenter + xlength; iphi++)
        {
          int index = (ieta - etacenter + ylength) * inputDimx + iphi - phicenter + xlength;
          if (index < 0 || index >= vectorSize) {
            std::cerr << "Mother index out of bounds: " << index << std::endl;
            std::cerr << "ieta: " << ieta << ", iphi: " << iphi << std::endl;
            continue; 
          }
          if(ieta < 0 || ieta >=96){
            input.at(index) = 0;
            continue;
          }
          int mappediphi = iphi;

          if (mappediphi < 0)
          {
            mappediphi += 256;
          }
          if (mappediphi > 255)
          {
            mappediphi -= 256;
          }
          unsigned int towerinfokey = TowerInfoDefs::encode_emcal(ieta, mappediphi);

          TowerInfo *towerinfo = offlinetowers->get_tower_at_key(towerinfokey);
          if (!towerinfo)
          {
            std::cout << "No towerinfo for tower key " << towerinfokey << std::endl;
            std::cout << "ieta: " << ieta << " iphi: " << mappediphi << std::endl;
            continue;
          }
          input.at(index) = towerinfo->get_energy();
          
          g_forfit->AddPoint(index/7,index%7,towerinfo->get_energy());
          if (towerinfo->get_energy() > 0.07) h_cluster->Fill(index/7,index%7,towerinfo->get_energy());
          if (towerinfo->get_energy() > emax1) 
          {
            emax2 = emax1;
            emax1 = towerinfo->get_energy();
          }
        }
      }
      nearbytowers_mother.push_back(input);
      
      TLorentzVector fitphoton1;
      TLorentzVector fitphoton2;
      TLorentzVector fitdiphoton;
      // Now try fitting
      int numBins = 7;
      if (clE > 7 && g_forfit->GetN() > 0) 
      {
        TH1D * cx = h_cluster->ProjectionX();
        TH1D * cy = h_cluster->ProjectionY();
        Double_t * xeta = new Double_t[numBins];
        Double_t * xphi = new Double_t[numBins];
        Double_t * yeta = new Double_t[numBins];
        Double_t * yphi = new Double_t[numBins];
        for (int i = 0; i < numBins; i++) {
          xeta[i] = cx->GetBinCenter(i);
          xphi[i] = cy->GetBinCenter(i);
          yeta[i] = cx->GetBinContent(i);
          yphi[i] = cy->GetBinContent(i);
        }
        double medianx = TMath::Median(numBins, xeta, yeta);
        double mediany = TMath::Median(numBins, xphi, yphi);
        float rmsx = h_cluster->GetStdDev(1);
        float rmsy = h_cluster->GetStdDev(2);
        float rms = TMath::Sqrt(TMath::Power(rmsx,2)+TMath::Power(rmsy,2));
        float meanx = h_cluster->GetMean(1);
        float meany = h_cluster->GetMean(2);
        int skewx = -1 + 2*(meanx - medianx > 0);
        int skewy = -1 + 2*(meany - mediany > 0);
        
        float corr = h_cluster->GetCorrelationFactor();
        int corrsign = corr/abs(corr);
        float eangle = TMath::ATan2(rmsy,corrsign*rmsx) + M_PI;

        float etacm = eta_cluster_mother_center[nClusters_mother];
        float phicm = phi_cluster_mother_center[nClusters_mother];
        float etacm_ff = etacm - (int)etacm + ( (etacm - (int)etacm > 0.5) ? 2 :  3 );
        float phicm_ff = phicm - (int)phicm + ( (phicm - (int)phicm > 0.5) ? 2 :  3 );

        float eta1_ff = etacm_ff + rms*TMath::Cos(eangle);
        float eta2_ff = etacm_ff - rms*TMath::Cos(eangle);
        float phi1_ff = phicm_ff + rms*TMath::Sin(eangle);
        float phi2_ff = phicm_ff - rms*TMath::Sin(eangle);
 
        float loweta1 =  (eta1_ff - 3 < 0) ? 0 : eta1_ff - 3;
        float loweta2 =  (eta2_ff - 3 < 0) ? 0 : eta2_ff - 3;
        float higheta1 = (eta1_ff + 3 > 7) ? 7 : eta1_ff + 3;
        float higheta2 = (eta2_ff + 3 > 7) ? 7 : eta2_ff + 3;
        float lowphi1 =  (phi1_ff - 3 < 0) ? 0 : phi1_ff - 3;
        float lowphi2 =  (phi2_ff - 3 < 0) ? 0 : phi2_ff - 3;
        float highphi1 = (phi1_ff + 3 > 7) ? 7 : phi1_ff + 3;
        float highphi2 = (phi2_ff + 3 > 7) ? 7 : phi2_ff + 3;

        float dr_ff = 2*rms;
        float dr_fa = dr_ff * 0.024;
        if (dr_fa < 0.004) dr_fa = 0.005; 
        float alpha_ff = 2*TMath::ACos(0.004/dr_fa)/M_PI;
        float ehigh_ff = clE*(1+alpha_ff)/2.0;
        float elow_ff = clE*(1-alpha_ff)/2.0;

        float skewangle = TMath::ATan2(skewx,skewy) + M_PI;
        int skewquad = (int) (skewangle / (M_PI/2.0)) + 1;
        int quad1 = (int) (eangle / (M_PI/2.0)) + 1;
        int quad2 = (quad1 + 2) % 4;

        if (skewquad == quad1)      { emax1 = elow_ff; emax2 = ehigh_ff; }
        else if (skewquad == quad2) { emax2 = elow_ff; emax1 = ehigh_ff; }
        if (emax1 < 0.1*clE) emax1 = 0.11*clE;
        if (emax1 > 0.9*clE) emax1 = 0.89*clE;
        if (emax2 < 0.1*clE) emax2 = 0.11*clE;
        if (emax2 > 0.9*clE) emax2 = 0.89*clE;
        
        TF2 * func = new TF2("func", [=, this](double *x, double *par) {
            Double_t r1 = TMath::Sqrt(TMath::Power(x[0]-par[1],2)+TMath::Power(x[1]-par[2],2));
            Double_t r2 = TMath::Sqrt(TMath::Power(x[0]-par[4],2)+TMath::Power(x[1]-par[5],2));
            Double_t value1;
            Double_t value2;
            if (r1 > 2.5) value1 = 0.0;
            else value1 = par[0]*profilehist->Interpolate(r1);
            if (r2 > 2.5) value2 = 0.0;
            else value2 = par[3]*profilehist->Interpolate(r2);
            return value1 + value2;
            }, -0.5,6.5,0.5,6.5, 6);
        
        func->SetParameters(emax1,eta1_ff,phi1_ff,emax2,eta2_ff,phi2_ff);
        func->SetParLimits(0,0.1*clE,0.9*clE);
        func->SetParLimits(1,loweta1,higheta1);
        func->SetParLimits(2,lowphi1,highphi1);
        func->SetParLimits(3,0.1*clE,0.9*clE);
        func->SetParLimits(4,loweta2,higheta2);
        func->SetParLimits(5,lowphi2,highphi2);

        // The star of the show
        fitres[nClusters_mother] = h_cluster->Fit("func","RQ0");

        // calculate cluster kinematics
        float fit_e1 = func->GetParameter(0);
        float fit_eta1 = clEta + (func->GetParameter(1)-etacm_ff) * 0.024;
        float fit_phi1 = clPhi + (func->GetParameter(2)-phicm_ff) * M_PI/128.0;
        float fit_e2 = func->GetParameter(3);
        float fit_eta2 = clEta + (func->GetParameter(4)-etacm_ff) * 0.024;
        float fit_phi2 = clPhi + (func->GetParameter(5)-phicm_ff) * M_PI/128.0;
        //float alpha = abs(fit_e1-fit_e2)/(fit_e1+fit_e2);
        //if (alpha > 0.9) continue;

        fitphoton1.SetPtEtaPhiE(fit_e1/TMath::CosH(fit_eta1),fit_eta1,fit_phi1,fit_e1);
        fitphoton2.SetPtEtaPhiE(fit_e2/TMath::CosH(fit_eta2),fit_eta2,fit_phi2,fit_e2);
        fitdiphoton = fitphoton1 + fitphoton2;

        new ((*photon1_4mom_fit)[nClusters_mother]) TLorentzVector(fitphoton1);
        new ((*photon2_4mom_fit)[nClusters_mother]) TLorentzVector(fitphoton2);
        new ((*diphoton_4mom_fit)[nClusters_mother]) TLorentzVector(fitdiphoton);
        fit_chisq[nClusters_mother] = func->GetChisquare();

        delete func;
        delete cx;
        delete cy;
        delete[] xeta;
        delete[] yeta;
        delete[] xphi;
        delete[] yphi;
      }
      else {
        fitphoton1.SetPtEtaPhiE(0,0,0,0);
        fitphoton2.SetPtEtaPhiE(0,0,0,0);
        fitdiphoton.SetPtEtaPhiE(0,0,0,0);
        new ((*photon1_4mom_fit)[nClusters_mother]) TLorentzVector(fitphoton1);
        new ((*photon2_4mom_fit)[nClusters_mother]) TLorentzVector(fitphoton2);
        new ((*diphoton_4mom_fit)[nClusters_mother]) TLorentzVector(fitdiphoton);
        fit_chisq[nClusters_mother] = -1;
      }

      delete h_cluster;
      delete g_forfit;

      nClusters_mother++;
    }
  
    


    //Subcluster loop
    nClusters=0;
    nPairs=0;
    int _hiter1=0;
    short photon1_idx = 0;

    nClusters=0;
    for(hiter = begin_end.first; hiter != begin_end.second; ++hiter)
    {
      _hiter1++;
      RawCluster* cluster = hiter->second;
      CLHEP::Hep3Vector vertex;
      if(isMC && UseClusterTruthVertex) vertex.set(truth_vx,truth_vy,truth_vz);
      else if(UseClusterVertexZero) vertex.set(0,0,0);
      else vertex.set(vx,vy,vz);
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetEVec(*cluster, vertex);
      
      float clPt = E_vec_cluster.perp();
      float clEta = E_vec_cluster.pseudoRapidity();
      float clPhi = E_vec_cluster.phi();
      float clE = E_vec_cluster.mag();
      //if(clPt < _emcal_cluster_emincut) continue;
      photon_prob[nClusters] = cluster->get_prob();
      photon_prob_merged_cluster[nClusters] = cluster->get_merged_cluster_prob();

      TLorentzVector vphoton1;
      vphoton1.SetPtEtaPhiE(clPt,clEta,clPhi,clE);

      photon1_idx++;
      short photon2_idx = 0;
      int _hiter2=0;
      for(hiter2 = begin_end.first; hiter2 != begin_end.second; ++hiter2){
        _hiter2++;
        if(hiter == hiter2) continue;
        if(_hiter1 >= _hiter2) continue;

        RawCluster* cluster2 = hiter2->second;
        CLHEP::Hep3Vector E_vec_cluster2 = RawClusterUtility::GetEVec(*cluster2, vertex);
        float clPt2 = E_vec_cluster2.perp();
        float clEta2 = E_vec_cluster2.pseudoRapidity();
        float clPhi2 = E_vec_cluster2.phi();
        float clE2 = E_vec_cluster2.mag();

        if( abs(clEta2-clEta)<1e-4 && abs(clPhi2-clPhi)<1e-4) continue;
        //if (clPt2 < _emcal_cluster_emincut) continue;

        TLorentzVector vphoton2;
        vphoton2.SetPtEtaPhiE(clPt2,clEta2,clPhi2,clE2);

        TLorentzVector vdiphoton = vphoton1+ vphoton2;
        new ((*diphoton_4mom)[nPairs]) TLorentzVector(vdiphoton);

        if(vdiphoton.M() < diphotonmasslow || vdiphoton.M()>diphotonmasshigh) continue;
        diphoton_energyimbal[nPairs] = fabs(clE - clE2) / (clE + clE2);
        diphoton_dR[nPairs] = vphoton1.DeltaR(vphoton2);

        photon2_idx++;
        idx_photon1[nPairs] = photon1_idx -1;
        idx_photon2[nPairs] = photon1_idx -1 + photon2_idx;
        
        if (isMC) {
          short closest_i = -1;
          float closest_dR = 999;
          for (short i = 0; i < truthpar_n; i++) {
            auto truthdiphov = (TLorentzVector*) truth_diphoton_4mom->At(i);
            float dphi = vdiphoton.Phi() - truthdiphov->Phi();
            if (dphi > M_PI) dphi -= M_PI;
            float deta = vdiphoton.Eta() - truthdiphov->Eta();
            float dR = TMath::Sqrt(dphi*dphi + deta*deta);
            if (dR < 0.1 && dR < closest_dR) {
              closest_i = i;
            }
          }
          reco_matched_truth_idx[nPairs] = closest_i;
        }
        nPairs++;
      }
      new ((*photon_4mom)[nClusters]) TLorentzVector(vphoton1);

      std::vector<float> showershape = cluster->get_shower_shapes(0.07);
      eta_cluster_center[nClusters] = showershape.at(4);
      phi_cluster_center[nClusters] = showershape.at(5);

      int etacenter = std::floor(eta_cluster_center[nClusters] + 0.5);
      int phicenter = std::floor(phi_cluster_center[nClusters] + 0.5);

      std::vector<float> input;
      int vectorSize = inputDimx * inputDimy;
      input.resize(vectorSize, 0);

      int xlength = int((inputDimx - 1) / 2);
      int ylength = int((inputDimy - 1) / 2);
      for (int ieta = etacenter - ylength; ieta <= etacenter + ylength; ieta++)
      {
        for (int iphi = phicenter - xlength; iphi <= phicenter + xlength; iphi++)
        {
          int index = (ieta - etacenter + ylength) * inputDimx + iphi - phicenter + xlength;
          if (index < 0 || index >= vectorSize) {
            std::cerr << "Index out of bounds: " << index << std::endl;
            std::cerr << "ieta: " << ieta << ", iphi: " << iphi << std::endl;
            continue; 
          }
          if(ieta < 0 || ieta >=96){
            input.at(index) = 0;
            continue;
          }
          int mappediphi = iphi;

          if (mappediphi < 0)
          {
            mappediphi += 256;
          }
          if (mappediphi > 255)
          {
            mappediphi -= 256;
          }
          unsigned int towerinfokey = TowerInfoDefs::encode_emcal(ieta, mappediphi);

          TowerInfo *towerinfo = offlinetowers->get_tower_at_key(towerinfokey);
          if (!towerinfo)
          {
            std::cout << "No towerinfo for tower key " << towerinfokey << std::endl;
            std::cout << "ieta: " << ieta << " iphi: " << mappediphi << std::endl;
            continue;
          }

          input.at(index) = towerinfo->get_energy();
        }
      }
      nearbytowers.push_back(input);
      nClusters++;
    }
  }
    //if(!isMC && nPairs==0) return Fun4AllReturnCodes::ABORTEVENT;
  
  if(doJets){

    //RawClusterContainer *clustersEM = static_cast<RawClusterContainer*>(findNode::getClass<RawClusterContainer>(topNode, "CEMC_CLUSTERINFO"));
    JetContainer* Jets = static_cast<JetContainer*>(findNode::getClass<JetContainer>(topNode, "AntiKt_unsubtracted_r04"));
    if(!Jets)
    {
      std::cout <<"JetNconstituents::process_event - Error can not find jet node " << "AntiKt_Tower_r04" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    TowerInfoContainer *towersEM3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
    TowerInfoContainer *towersIH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
    TowerInfoContainer *towersOH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT"); 
    if(!towersEM3 || !towersIH3 || !towersOH3)
    {
      std::cout <<"JetNconstituents::process_event - Error can not find tower node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    // get tower geometry
    RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    if(!geomIH || !geomOH)
    {
      std::cout <<"JetNconstituents::process_event - Error can not find tower geometry node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    nJets=0;
    for (auto jet : *Jets)
    {    
      if(jet->get_e() < _jet_minetcut) continue;
      if(!(fabs(jet->get_eta()) < _jet_etacut)) continue;

      jet_e[nJets] = jet->get_e();
      jet_pt[nJets] = jet->get_pt();
      jet_eta[nJets] = jet->get_eta();
      jet_phi[nJets] = jet->get_phi();

      float m_jet_total_eT = 0;
      for (auto comp: jet->get_comp_vec())
      {
        unsigned int channel = comp.second;
        TowerInfo *tower;

        if(comp.first == 26 || comp.first == 30)
        { // IHcal 

          tower = towersIH3->get_tower_at_channel(channel);
          if(!tower || !geomIH){ continue; }
          if(!tower->get_isGood()) continue;

          unsigned int calokey = towersIH3->encode_key(channel);
          int ieta = towersIH3->getTowerEtaBin(calokey);
          int iphi = towersIH3->getTowerPhiBin(calokey);
          const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
          //float tower_phi = tower_geomIH->get_tower_geometry(key)->get_phi();
          float tower_eta = geomIH->get_tower_geometry(key)->get_eta();
          float tower_eT = tower->get_energy()/cosh(tower_eta);

          m_jetfrac_ihcal[nJets]+= tower_eT;
          m_jet_total_eT += tower_eT;
        }
        else if(comp.first == 27 || comp.first == 31)
        { // OHcal 

          tower = towersOH3->get_tower_at_channel(channel);

          if(!tower || !geomOH){ continue; }
          if(!tower->get_isGood()) continue;

          unsigned int calokey = towersOH3->encode_key(channel);
          int ieta = towersOH3->getTowerEtaBin(calokey);
          int iphi = towersOH3->getTowerPhiBin(calokey);
          const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);
          //float tower_phi = tower_geomOH->get_tower_geometry(key)->get_phi();
          float tower_eta = geomOH->get_tower_geometry(key)->get_eta();
          float tower_eT = tower->get_energy()/cosh(tower_eta);
          m_jetfrac_ohcal[nJets] += tower_eT;
          m_jet_total_eT += tower_eT;
        }
        else if(comp.first == 28 || comp.first == 29)
        { // EMCal (retowered)

          tower = towersEM3->get_tower_at_channel(channel);

          if(!tower || !geomIH){ continue; }
          if(!tower->get_isGood()) continue;

          unsigned int calokey = towersEM3->encode_key(channel);
          int ieta = towersEM3->getTowerEtaBin(calokey);
          int iphi = towersEM3->getTowerPhiBin(calokey);
          const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
          //float tower_phi = tower_geomIH->get_tower_geometry(key)->get_phi();
          float tower_eta = geomIH->get_tower_geometry(key)->get_eta();
          float tower_eT = tower->get_energy()/cosh(tower_eta);

          m_jetfrac_emcal[nJets] += tower_eT;
          m_jet_total_eT += tower_eT;
        }
      }

      m_jetfrac_ihcal[nJets] /= m_jet_total_eT;
      m_jetfrac_ohcal[nJets] /= m_jet_total_eT;
      m_jetfrac_emcal[nJets] /= m_jet_total_eT;
      nJets++;
    }    
    if(isMC){
      JetContainer* TruthJets = findNode::getClass<JetContainer>(topNode, "AntiKt_Truth_r04");
      if(!TruthJets)
      {
        std::cout <<"Error can not find jet truth node " << "AntiKt_Truth_r04" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
      nTruthJets=0;
      for (auto truthjet : *TruthJets)
      {    
        if(truthjet->get_e() < _jet_minetcut) continue;
        if(!(fabs(truthjet->get_eta()) < _jet_etacut)) continue;

        truth_jet_e[nTruthJets] = truthjet->get_e();
        truth_jet_pt[nTruthJets] = truthjet->get_pt();
        truth_jet_eta[nTruthJets] = truthjet->get_eta();
        truth_jet_phi[nTruthJets] = truthjet->get_phi();
        nTruthJets++;
      }
    }
  }  

  if(nPairsMother==1 && nClusters_mother==1) 
  {
    std::cout << "something is very wrong..." << std::endl;
  }


  if(doClusters && doJets){
    if(nClusters ==0 || nJets==0) return Fun4AllReturnCodes::ABORTEVENT;
  }

  //if (doClusters && isMC) {
  //  if (nPairs == 0 && truthpar_n == 0) return Fun4AllReturnCodes::ABORTEVENT;
  //}

  towerntuple->Fill();


  m_emcal_total_energy = 0;
  m_mbd_energy_south = 0;
  m_mbd_energy_north = 0;
  if(doClusters){
    photon_4mom->Clear();
    photon_4mom_mother->Clear();
    diphoton_4mom->Clear();
    diphoton_4mom_mother->Clear();
    nearbytowers.clear();
    nearbytowers_mother.clear();
    photon1_4mom_fit->Clear();
    photon2_4mom_fit->Clear();
    diphoton_4mom_fit->Clear();
    if (isMC) {
      truth_photon1_4mom->Clear();
      truth_photon2_4mom->Clear();
      truth_diphoton_4mom->Clear();
    }

  }

  return Fun4AllReturnCodes::EVENT_OK;
}
  

int CaloAna::ProcessGlobalEventInfo(PHCompositeNode* topNode){

  vz=-999;
  vx=-999;
  vy=-999;
    

  if(isMC){
    if(useEmulator){
      triggeranalyzer->decodeTriggers(topNode);
      ScaledTriggerBit[20] =  triggeranalyzer->didTriggerFire("Jet 6 GeV");
      ScaledTriggerBit[21] =  triggeranalyzer->didTriggerFire("Jet 8 GeV");
      ScaledTriggerBit[22] =  triggeranalyzer->didTriggerFire("Jet 10 GeV");
      ScaledTriggerBit[23] =  triggeranalyzer->didTriggerFire("Jet 12 GeV");
      ScaledTriggerBit[28] =  triggeranalyzer->didTriggerFire("Photon 2 GeV");
      ScaledTriggerBit[29] =  triggeranalyzer->didTriggerFire("Photon 3 GeV");
      ScaledTriggerBit[30] =  triggeranalyzer->didTriggerFire("Photon 4 GeV");
      ScaledTriggerBit[31] =  triggeranalyzer->didTriggerFire("Photon 5 GeV");
    }
    MbdVertexMap *mbdvtxmap = findNode::getClass<MbdVertexMap>(topNode,"MbdVertexMap");
    bool isglbvtx=true;
    if(!mbdvtxmap || mbdvtxmap->empty()){ 
      //std::cout << "Empty mbdmap or mbdvtx node" << std::endl;
      isglbvtx=false;
    }
    if(isglbvtx){
      MbdVertex *bvertex= nullptr;
      if (mbdvtxmap && m_UseAltZVertex == 1)
      {
        for (MbdVertexMap::ConstIter mbditer= mbdvtxmap->begin(); mbditer != mbdvtxmap->end(); ++mbditer)
        {
          bvertex = mbditer->second;
        }
        if(!bvertex){std::cout << "could not find globalvtxmap iter :: set vtx to (-999,-999,-999)" << std::endl;}
        else if(bvertex){
          vz = bvertex->get_z();
          vy = bvertex->get_y();
          vx = bvertex->get_x();
          count_mbdvtx++;
        }
      }
    }

    PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (truthinfo)
    {
      PHG4VtxPoint *gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
      truth_vz = gvertex->get_z();
      truth_vx = gvertex->get_x();
      truth_vy = gvertex->get_y();
    }
  }
  else if(!isMC){
    Gl1Packet *gl1_packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (gl1_packet)
    {
      uint64_t gl1_scaledtriggervector = gl1_packet->lValue(0, "ScaledVector");
      uint64_t gl1_livetriggervector = gl1_packet->lValue(0, "TriggerVector");

      for (int i = 0; i < 64; i++)
      {
        ScaledTriggerBit[i] = ((gl1_scaledtriggervector & 0x1U) == 0x1U);
        LiveTriggerBit[i] = ((gl1_livetriggervector & 0x1U) == 0x1U);
        count_raw[i] = gl1_packet->lValue(i, 0);
        count_live[i] = gl1_packet->lValue(i, 1);    
        count_scaled[i] = gl1_packet->lValue(i, 2);
        gl1_scaledtriggervector = (gl1_scaledtriggervector >> 1U) & 0xffffffffU;
        gl1_livetriggervector = (gl1_livetriggervector >> 1U) & 0xffffffffU;
      }
      if(LiveTriggerBit[10] == true) mbdlivecount++;
    }
    GlobalVertexMap *globalvtxmap = findNode::getClass<GlobalVertexMap>(topNode,"GlobalVertexMap");
    bool isglbvtx=true;
    if(!globalvtxmap) return Fun4AllReturnCodes::ABORTEVENT;
    if(globalvtxmap->empty()){ 
      isglbvtx=false;
    //  return Fun4AllReturnCodes::ABORTEVENT;
    }
    if(isglbvtx){
      GlobalVertex *bvertex= nullptr;
      if (globalvtxmap && m_UseAltZVertex == 1)
      {
        for (GlobalVertexMap::ConstIter globaliter= globalvtxmap->begin(); globaliter != globalvtxmap->end(); ++globaliter)
        {
          bvertex = globaliter->second;
        }
        if(!bvertex){std::cout << "could not find globalvtxmap iter :: set vtx to (-999,-999,-999)" << std::endl;}
        else if(bvertex){
          vz = bvertex->get_z();
          vy = bvertex->get_y();
          vx = bvertex->get_x();
          count_mbdvtx++;
        }
      }
    }
  }

  //Event info
  if(isMC){
    if(isPythia){
      PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
      if(!genevtmap){std::cout << "no PHHepMCGenEventMap" << std::endl; return Fun4AllReturnCodes::ABORTEVENT;}
      PHHepMCGenEvent *genevt = nullptr;

      //for (PHHepMCGenEventMap::ReverseIter iter = genevtmap->rbegin(); iter != genevtmap->rend(); ++iter)
      //{   
        genevt = genevtmap->rbegin()->second;
        //genevt = iter->second;
        if(!genevt) {cout<<"ERROR: no PHHepMCGenEvent!" << endl; return Fun4AllReturnCodes::ABORTEVENT;}
      //}   

      HepMC::GenEvent *event = genevt->getEvent();
      if (!event) { cout << PHWHERE << "ERROR: no HepMC::GenEvent!" << endl; return Fun4AllReturnCodes::ABORTEVENT;}
      processId = event->signal_process_id();
    }
    if(isHI){
      EventHeaderv1 *event_header = findNode::getClass<EventHeaderv1>(topNode, "EventHeader" );
      if ( event_header ) {
        npart = event_header->get_intval("npart");
        ncoll = event_header->get_intval("ncoll");
        bimp = event_header->get_floatval("bimp");
      } 
      CentralityInfo* cent_node = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
      if (!cent_node)
      {
        std::cout
          << "Error can not find centrality node "
          << std::endl;
        exit(-1);
      }   
      Centrality =  cent_node->get_centile(CentralityInfo::PROP::mbd_NS);
      Cent_impactparam =  cent_node->get_quantity(CentralityInfo::PROP::bimp);
    }
  }
  else if(!isMC){
    if(isHI){
      CentralityInfo *centrality = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
      if (!centrality)
      {
        std::cout << "no centrality node " << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }

      MinimumBiasInfo *minimumbiasinfo = findNode::getClass<MinimumBiasInfo>(topNode, "MinimumBiasInfo");
      if (!minimumbiasinfo)
      {
        std::cout << "no minimumbias node " << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }

      float centile = (centrality->has_centile(CentralityInfo::PROP::mbd_NS) ? centrality->get_centile(CentralityInfo::PROP::mbd_NS) : -999.99);
      centbin = centile*100;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloAna::ProcessFillTruthParticle(PHCompositeNode* topNode, float truthvz){
  PHG4TruthInfoContainer* truthinfo = findNode::getClass <PHG4TruthInfoContainer> (topNode, "G4TruthInfo");
  if(!truthinfo){std::cout << "no truth info node... just skip the whole part.." << std::endl; return;}
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  
  truth_photon1_4mom->Clear();
  truth_photon2_4mom->Clear();
  truth_diphoton_4mom->Clear(); 
    
  float etamin = GetShiftedEta(truthvz,-2);
  float etamax = GetShiftedEta(truthvz,2);
  truthpar_n=0;
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    PHG4Particle* g4particle = iter->second;
    if(isSingleGun){
      if(g4particle->get_pid() != PDGPID || g4particle->get_parent_id()!=0 || !truthinfo->is_primary(g4particle)){
        std::cout << "something is not right... " << g4particle->get_pid() << "/" << PDGPID << " " << g4particle->get_parent_id() << " " << truthinfo->is_primary(g4particle) << std::endl;
      }
    }
    if (!truthinfo->is_primary(g4particle)) continue;
    if(isHI){
      if(truthinfo->isEmbeded(g4particle->get_track_id())<=0) continue;
    }

    if(g4particle->get_parent_id()!=0) continue;
    if(isPythia){
      //if(g4particle->get_pid() != 111 && g4particle->get_pid() != 221 && g4particle->get_pid() !=22 ) continue;
      if(g4particle->get_pid() != 111 && g4particle->get_pid() != 221 ) continue;
    }
    else if(!isPythia){
      if(g4particle->get_pid() != PDGPID || g4particle->get_parent_id()!=0) continue;
    }

    TLorentzVector t;
    t.SetPxPyPzE(g4particle->get_px (), g4particle->get_py (), g4particle->get_pz (), g4particle->get_e ());
    new ((*truth_diphoton_4mom)[truthpar_n]) TLorentzVector(t);

    if(t.Eta() > etamax || t.Eta() < etamin) continue;
    // take only particles from primary Pythia event                                                            
    truth_id[truthpar_n] = g4particle->get_pid();
    isconverted[truthpar_n]= FindConversion(truthinfo,g4particle->get_track_id(),t.E());
    truth_track_id[truthpar_n]= g4particle->get_track_id();
    found_decay1[truthpar_n] = false;
    found_decay2[truthpar_n] = false;
    isconverted1[truthpar_n]= false;
    isconverted2[truthpar_n]= false;
    new ((*truth_photon1_4mom)[truthpar_n]) TLorentzVector(0, 0, 0, 0);
    new ((*truth_photon2_4mom)[truthpar_n]) TLorentzVector(0, 0, 0, 0);
    truthpar_n++;
  }

  if(!isPythia && PDGPID == 22) return;
  PHG4TruthInfoContainer::Range truthRangeDecay1 = truthinfo->GetSecondaryParticleRange();
  PHG4TruthInfoContainer::ConstIterator truthIterDecay1;
  for(truthIterDecay1 = truthRangeDecay1.first; truthIterDecay1 != truthRangeDecay1.second; truthIterDecay1++)
  {
    PHG4Particle *decay = truthIterDecay1 -> second;
    int truthpid = decay -> get_pid();
    int parentid = decay -> get_parent_id();

    if(truthpid == 22){
      for(int it = 0; it<truthpar_n; it++){
        if(parentid != truth_track_id[it]){
          continue;
        }
        TLorentzVector tpho;
        tpho.SetPxPyPzE(decay -> get_px(), decay -> get_py(), decay -> get_pz(), decay -> get_e());

        if(!found_decay1[it]){
          new ((*truth_photon1_4mom)[it]) TLorentzVector(tpho);
          found_decay1[it] = true;
          isconverted1[it]= FindConversion(truthinfo,truth_track_id[it],tpho.E());
        }
        else if(!found_decay2[it]){
          new ((*truth_photon2_4mom)[it]) TLorentzVector(tpho);
          found_decay2[it] = true;
          isconverted2[it]= FindConversion(truthinfo,truth_track_id[it],tpho.E());
        }
      }
    }
  }
}

bool CaloAna::FindConversion(PHG4TruthInfoContainer* truthinfo, int trackid, float energy){
  PHG4Shower *shower = truthinfo->GetShower(trackid);
  if(!shower){
    return false;
  }
  bool foundconversion=false;
  auto g4particle_ids = shower->g4particle_ids();
  for (auto g4particle_id : g4particle_ids)
  {
    PHG4Particle *g4particle = truthinfo->GetParticle(g4particle_id);
    if (!g4particle) continue;
    int vertexid = g4particle->get_vtx_id();
    PHG4VtxPoint *vtxp = truthinfo->GetVtx(vertexid);
    if (!vtxp) continue;
    float vertexr = sqrt(vtxp->get_x() * vtxp->get_x() + vtxp->get_y() * vtxp->get_y());
    if (vertexr < 93)
    {
      float momentum = sqrt(g4particle->get_px() * g4particle->get_px() + g4particle->get_py() * g4particle->get_py() + g4particle->get_pz() * g4particle->get_pz());
      if (momentum > 0.3 * energy)
      {
        int g4particlepid = g4particle->get_pid();
        if (abs(g4particlepid) == 11) foundconversion=true;
      }
    }
  }
  return foundconversion;
}
  
Double_t CaloAna::GetShiftedEta(float _vz, float _eta)
{
  double radius = 93;
  double theta = 2*atan(exp(-_eta));
  double z = radius / tan(theta);
  double zshifted = z - _vz;
  double thetashifted = atan2(radius,zshifted);
  double etashifted = -log(tan(thetashifted/2.0));
  return etashifted;
}


int CaloAna::End(PHCompositeNode* /*topNode*/)
{
  if(!isMC){
    std::cout << "fill MBD live counts : " << mbdlivecount << std::endl;
    histlumicount->SetBinContent(1,mbdlivecount);
  }
  std::cout << "end..! " << std::endl;
  std::cout << "Events with MBD N&S>=1 live bit raised : " << mbdlivecount << " / events with a valid MBD vertex : " << count_mbdvtx << " / events with MBD N&S hit>=1 " << count_mbdhit << " out of total events : " << count << std::endl;
  if(outfile){
    outfile->cd();
    towerntuple->Write();
    if(!isMC) histlumicount->Write();
    outfile->Write();
    outfile->Close();
    delete outfile;
    outfile=nullptr;
  }
  return 0;
}
