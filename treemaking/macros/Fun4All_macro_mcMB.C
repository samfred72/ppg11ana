#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <phool/recoConsts.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>


#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>

#include <Calo_Calib.C>
#include <G4_Global.C>
#include <GlobalVariables.C>
#include <mbd/MbdReco.h>
#include <globalvertex/GlobalVertexReco.h>
//#include <calotowerbuilder/CaloTowerBuilder.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloWaveformProcessing.h>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <clusteriso/ClusterIso.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>

// #include <runtowerinfo/RunTowerInfo.h>
#include <caloana/CaloAna.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtContainerV1.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertexReco.h>
#include <mbd/MbdReco.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4centrality/PHG4CentralityReco.h>

#include <centrality/CentralityReco.h>
#include <calotrigger/MinimumBiasClassifier.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libcaloana.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4vertex.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libclusteriso.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libFROG.so)


#endif

void Fun4All_macro_mcMB(const char * inqueue="queue_run28_v00001.list", bool isMC=true, bool isEmbed=false, float mbdzvtxcut = -999)
{
    string infile = Form("/sphenix/user/samfred/run25/ppg11/treemaking/macros/filelists/pythia_MB/%s",inqueue);
    cout << "Infile: " << infile << endl;
    vector<string> infile_cluster;
    vector<string> infile_global;
    vector<string> infile_truth;
    vector<string> infile_mbd;
    
    ifstream file(infile); 
    std::string line;
    while (std::getline(file, line)) {
       std::istringstream iss(line);
       string dst1, dst2, dst3, dst4;
       if (iss >> dst1 >> dst2 >> dst3 >> dst4) {
           infile_cluster.push_back(dst1);
           infile_global.push_back(dst2);
           infile_truth.push_back(dst3);
           infile_mbd.push_back(dst4);
       }
    }
    
    Fun4AllServer *se = Fun4AllServer::instance();
    int verbosity = 0;

    se->Verbosity(verbosity);
    recoConsts *rc = recoConsts::instance();

    pair<int, int> runseg = Fun4AllUtils::GetRunSegment(infile_cluster.at(0));
    int runnumber = runseg.first;
    int segment = runseg.second;
    
    string outdir = Form("/sphenix/tg/tg01/jets/samfred/run25/pythia_MB");
    void * dirf = gSystem->OpenDirectory(outdir.c_str());
    if(dirf) gSystem->FreeDirectory(dirf);
    else {gSystem->mkdir(outdir.c_str(), kTRUE);}
    
    //const char *outfile = Form("%s/outtree_%s",outdir.c_str(),inqueue);
    const char *outfile = Form("test_%s.root",inqueue);

    //===============
    // conditions DB flags
    //===============

    // global tag
    rc->set_StringFlag("CDB_GLOBALTAG","MDC2"); 
    rc->set_uint64Flag("TIMESTAMP",runnumber);

    //===============
    // conditions DB flags
    //===============
    
    Process_Calo_Calib();

    Fun4AllInputManager *incluster = new Fun4AllDstInputManager("DST_CLUSTER");
    for (int i = 0; i < infile_cluster.size(); i++) {
      incluster->AddFile(infile_cluster.at(i));
      cout << infile_cluster.at(i) << endl;
    }
    se->registerInputManager(incluster);

    Fun4AllInputManager *inglobal = new Fun4AllDstInputManager("DST_GLOBAL");
    for (int i = 0; i < infile_global.size(); i++) {
      inglobal->AddFile(infile_global.at(i));
    }
    se->registerInputManager(inglobal);

    Fun4AllInputManager *inmbd = new Fun4AllDstInputManager("DST_MBD");
    for (int i = 0; i < infile_mbd.size(); i++) {
      inmbd->AddFile(infile_mbd.at(i));
    }
    se->registerInputManager(inmbd);

    Fun4AllInputManager *intruth = new Fun4AllDstInputManager("DST_TRUTH");
    for (int i = 0; i < infile_truth.size(); i++) {
      intruth->AddFile(infile_truth.at(i));
    }
    se->registerInputManager(intruth);

    CaloAna *ca = new CaloAna("caloana",outfile,1,0);
    ca->SetRunNumber(runnumber);
    ca->SetMbdZVtxCut(mbdzvtxcut);
    ca->SetMCFlags(1,0);
    ca->SetUseOfClusters(0);//should be 1
    se->registerSubsystem(ca);

    std::cout << "now run..." << std::endl;
    se->run(10000);
    se->End();
    std::cout << "ok done.. " << std::endl;

    delete se;
    
}
