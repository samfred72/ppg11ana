#include <iostream>
#include "EventFilter.h"
using namespace std;

  EventFilter::EventFilter(const std::string& name, const char* outname, bool _isMC, bool _isHI)
: SubsysReco(name)
  , detector("HCALIN")
{
  outfilename = Form("%s",outname);
  isMC = _isMC;
  isHI = _isHI;
}

EventFilter::~EventFilter()
{
    //delete hm;
    delete towerntuple;
}

int EventFilter::Init(PHCompositeNode*)
{ 
  return 0;
}


int EventFilter::process_event(PHCompositeNode* topNode)
{
    if(!topNode){
      std::cerr << "EventFilter::Init - topnode PHCompositeNode not valid! return -1" << std::endl;
      return -1; 
    }
    if(count % 1000 == 0) std::cout << "event : " << count  << std::endl;
    count++;
    process_towers(topNode);
    return Fun4AllReturnCodes::EVENT_OK;
}

int EventFilter::process_towers(PHCompositeNode* topNode)
{
  ProcessGlobalEventInfo(topNode);
  ProcessFillTruthParticle(topNode);
  if(truthpar_n ==0) return Fun4AllReturnCodes::ABORTEVENT;

  return Fun4AllReturnCodes::EVENT_OK;
}

int EventFilter::ProcessGlobalEventInfo(PHCompositeNode* topNode){

  //Trigger
  truth_vz=-999;
  truth_vx=-999;
  truth_vy=-999;

  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (truthinfo)
  {
    PHG4VtxPoint *gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
    truth_vz = gvertex->get_z();
    truth_vx = gvertex->get_x();
    truth_vy = gvertex->get_y();
  }


  //std::cout << "----------------------------" << std::endl;
  //std::cout << "filtering vz : " << truth_vz << std::endl;
  if(fabs(truth_vz) > 10 ) return Fun4AllReturnCodes::ABORTEVENT;
  std::cout << std::endl;
  std::cout << "******* pass vz : " << truth_vz << " ********** " << std::endl;
  
  //Event info
  return Fun4AllReturnCodes::EVENT_OK;
}

void EventFilter::ProcessFillTruthParticle(PHCompositeNode* topNode){
    PHG4TruthInfoContainer* truthinfo = findNode::getClass <PHG4TruthInfoContainer> (topNode, "G4TruthInfo");
    if(!truthinfo){std::cout << "no truth info node... just skip the whole part.." << std::endl; return;}
    PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
    truthpar_n=0;
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
    {

      PHG4Particle* g4particle = iter->second;
      if(g4particle->get_pid() != PDGPID || g4particle->get_parent_id()!=0 || !truthinfo->is_primary(g4particle)){
        std::cout << "something is not right... " << g4particle->get_pid() << "/" << PDGPID << " " << g4particle->get_parent_id() << " " << truthinfo->is_primary(g4particle) << std::endl;
      }
      if (!truthinfo->is_primary(g4particle)) continue;

      if(isHI){
        if(truthinfo->isEmbeded(g4particle->get_track_id())<=0) continue;
      }

      if(g4particle->get_pid() != PDGPID || g4particle->get_parent_id()!=0) continue;

      TLorentzVector t;
      t.SetPxPyPzE (g4particle->get_px (), g4particle->get_py (), g4particle->get_pz (), g4particle->get_e ());

      // take only particles from primary Pythia event                                                            
      truth_energy[truthpar_n] = t.E();
      truth_pt[truthpar_n] = t.Pt();
      truth_eta[truthpar_n] = t.Eta();
      truth_phi[truthpar_n] = t.Phi();
      truth_id[truthpar_n] = g4particle->get_pid();

      if(fabs(truth_eta[truthpar_n]) > 0.9 ) continue;
      truthpar_n++;
    }
}

int EventFilter::End(PHCompositeNode* /*topNode*/)
{
  std::cout << "---- Filter end ----" << std::endl;
  return 0;
}
