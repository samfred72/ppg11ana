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

//#include <truthneutralmesonfinder/TruthNeutralMeson.h>
//#include <truthneutralmesonfinder/TruthNeutralMesonv1.h>
//#include <truthneutralmesonfinder/TruthNeutralMesonContainer.h>
//#include <truthneutralmesonfinder/TruthNeutralMesonBuilder.h>

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
//R__LOAD_LIBRARY(libtruthneutralmeson.so)

#endif

void Fun4All_macro(const char * inqueue="queue_run28_v00001.list", bool isMC=true, bool isEmbed=false, float mbdzvtxcut = -999)
{
    string infile = Form("/sphenix/user/samfred/run25/ppg11/data_treemaking/macros/filelists/queue/%s",inqueue);
    cout << "Infile: " << infile << endl;
    vector<string> infile_dst;
    
    ifstream file(infile); 
    std::string line;
    while (std::getline(file, line)) {
       std::istringstream iss(line);
       string dst;
       if (iss >> dst) {
           infile_dst.push_back(dst);
       }
    }
   
    Fun4AllServer *se = Fun4AllServer::instance();
    int verbosity = 0;

    se->Verbosity(verbosity);
    recoConsts *rc = recoConsts::instance();

    pair<int, int> runseg = Fun4AllUtils::GetRunSegment(infile_dst.at(0));
    int runnumber = runseg.first;
    int segment = runseg.second;
    
    string outdir = Form("/sphenix/tg/tg01/jets/samfred/run25/ppg11ana");
    void * dirf = gSystem->OpenDirectory(outdir.c_str());
    if(dirf) gSystem->FreeDirectory(dirf);
    else {gSystem->mkdir(outdir.c_str(), kTRUE);}
    
    const char *outfile = Form("%s/outtree_%s",outdir.c_str(),inqueue);
    //const char *outfile = Form("test_%s.root",inqueue);

    //===============
    // conditions DB flags
    //===============

    // global tag
    rc->set_StringFlag("CDB_GLOBALTAG","ProdA_2024"); 
    cout << runnumber << endl;
    rc->set_uint64Flag("TIMESTAMP",runnumber);
    
    std::unique_ptr<FlagHandler> flag = std::make_unique<FlagHandler>();
    se->registerSubsystem(flag.release());

    //===============
    // conditions DB flags
    //===============
    
    MbdReco *mbdreco = new MbdReco();
    se->registerSubsystem(mbdreco);
    GlobalVertexReco *gvertex = new GlobalVertexReco();
    se->registerSubsystem(gvertex);
    Process_Calo_Calib();
    
    Fun4AllInputManager *in = new Fun4AllDstInputManager("DST");
    for (int i = 0; i < infile_dst.size(); i++) {
      in->AddFile(infile_dst.at(i));
      cout << infile_dst.at(i) << endl;
    }
    se->registerInputManager(in);

    CaloAna *ca = new CaloAna("caloana",outfile,0,0);
    ca->SetRunNumber(runnumber);
    ca->SetMbdZVtxCut(mbdzvtxcut);
    ca->SetMCFlags(0,0);
    ca->SetUseOfClusters(1);//should be 1
    se->registerSubsystem(ca);

    cout << "Printing CDB info..." << endl;
    CDBInterface::instance()->Verbosity(Fun4AllBase::VERBOSITY_SOME);
    CDBInterface::instance()->Print();
    std::cout << "now run..." << std::endl;
    se->run();
    se->End();
    std::cout << "ok done.. " << std::endl;

    delete se;
    
}
