#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I/afs/rhic.bnl.gov/eic/include/ -I/afs/rhic.bnl.gov/eic/env/new/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent createJetTrees_aktR1.0Pt_v3.cxx -o createJetTrees_aktR1.0Pt_v3.exe

#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I/afs/rhic.bnl.gov/eic/include/ -I/afs/rhic.bnl.gov/eic/env/new/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent createJetTreesSmeared_aktR1.0_v3.cxx -o createJetTreesSmeared_aktR1.0_v3.exe

#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I/afs/rhic.bnl.gov/eic/include/ -I/afs/rhic.bnl.gov/eic/env/new/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent createJetTrees_partonRad_v3.cxx -o createJetTrees_partonRad_v3.exe

#echo $EICDIRECTORY

#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I$EICDIRECTORY/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent createJetTrees_testWTA.cxx -o createJetTrees_testWTA.exe

#g++ -Wall -ggdb -lLHAPDF `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I$EICDIRECTORY/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -L$EICDIRECTORY/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -leicsmear -lStEpSimuJetEvent createJetTrees_aktR1.0Pt_v3.cxx -o createJetTrees_aktR1.0Pt_v3.exe

g++ -Wall -ggdb -lLHAPDF `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I$EICDIRECTORY/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -L$EICDIRECTORY/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -leicsmear -lStEpSimuJetEvent createJetTrees_heraJetShape_v3.cxx -o createJetTrees_heraJetShape_v3.exe

#g++ -Wall -ggdb -lLHAPDF `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I$EICDIRECTORY/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -L$EICDIRECTORY/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -leicsmear -lStEpSimuJetEvent createJetTrees_angRadPt_v3.cxx -o createJetTrees_angRadPt_v3.exe

#g++ -Wall -ggdb -lLHAPDF `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I$EICDIRECTORY/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -L$EICDIRECTORY/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -leicsmear -lStEpSimuJetEvent createJetTreesSmeared_aktR1.0Pt_v3.cxx -o createJetTreesSmeared_aktR1.0Pt_v3.exe

#g++ -Wall -ggdb -lLHAPDF `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I$EICDIRECTORY/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -L$EICDIRECTORY/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -leicsmear -lStEpSimuJetEvent createJetTrees_angRadPt_v3.cxx -o createJetTrees_angRadPt_v3.exe

#g++ -Wall -ggdb -lLHAPDF `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I$EICDIRECTORY/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -L$EICDIRECTORY/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -leicsmear -lStEpSimuJetEvent createJetTrees_algoRadTestPt_v3.cxx -o createJetTrees_algoRadTestPt_v3.exe

#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I/afs/rhic.bnl.gov/eic/include/ -I/afs/rhic.bnl.gov/eic/env/new/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent createJetTrees_algoRadTestPt_v3.cxx -o createJetTrees_algoRadTestPt_v3.exe

#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I/afs/rhic.bnl.gov/eic/include/ -I/afs/rhic.bnl.gov/eic/env/new/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent createJetTrees_algoRadTest_v3.cxx -o createJetTrees_algoRadTest_v3.exe

#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I/afs/rhic.bnl.gov/eic/include/ -I/afs/rhic.bnl.gov/eic/env/new/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent createJetTrees_boostTest_v3.cxx -o createJetTrees_boostTest_v3.exe

#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I/afs/rhic.bnl.gov/eic/include/ -I/afs/rhic.bnl.gov/eic/env/new/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent createJetTrees_eektTest.cxx -o createJetTrees_eektTest.exe

#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I/afs/rhic.bnl.gov/eic/include/ -I/afs/rhic.bnl.gov/eic/env/new/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent testJetTree.cxx -o testJetTree.exe

#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I/afs/rhic.bnl.gov/eic/include/ -I/afs/rhic.bnl.gov/eic/env/new/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent createJetTrees_aktR1.0_v2.cxx -o createJetTrees_aktR1.0_v2.exe


#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I/afs/rhic.bnl.gov/eic/include/ -I/afs/rhic.bnl.gov/eic/env/new/include/ -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent createJetTrees_aktR1.0_v1.cxx -o createJetTrees_aktR1.0_v1.exe



#g++ -Wall -ggdb `root-config --cflags --libs` -I/eic/u/bpage/epJets/StEpSimuJetMaker -L/eic/u/bpage/epJets/StEpSimuJetMaker -lStEpSimuJetEvent analyzeJetTrees_v1.cxx -o analyzeJetTrees_v1.exe

#g++ -Wall -ggdb `root-config --cflags --libs` -I/eic/u/bpage/epJets/StEpSimuJetMaker -L/eic/u/bpage/epJets/StEpSimuJetMaker -lStEpSimuJetEvent analyzeJetTreesRemnant_v1.cxx -o analyzeJetTreesRemnant_v1.exe

#g++ -Wall -ggdb `root-config --cflags --libs` -I/eic/u/bpage/epJets/StEpSimuJetMaker -L/eic/u/bpage/epJets/StEpSimuJetMaker -lStEpSimuJetEvent analyzeJetTreesKTTest_v1.cxx -o analyzeJetTreesKTTest_v1.exe

#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I/afs/rhic.bnl.gov/eic/include/ -I/afs/rhic.bnl.gov/eic/env/new/include/ -I/eic/u/bpage/epJets/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent createJetTrees_v1.cxx -o createJetTrees_v1.exe

#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I/afs/rhic.bnl.gov/eic/include/ -I/afs/rhic.bnl.gov/eic/env/new/include/ -I/eic/u/bpage/epJets/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent createJetTrees_v1a.cxx -o createJetTrees_v1a.exe

#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I/afs/rhic.bnl.gov/eic/include/ -I/afs/rhic.bnl.gov/eic/env/new/include/ -I/eic/u/bpage/epJets/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent testHBFrame.cxx -o testHBFrame.exe

#g++ -Wall -ggdb `root-config --cflags --libs` `fastjet-config --cxxflags --libs --plugins` -I/afs/rhic.bnl.gov/eic/include/ -I/afs/rhic.bnl.gov/eic/env/new/include/ -I/eic/u/bpage/epJets/StEpSimuJetMaker -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/bpage/epJets/StEpSimuJetMaker -leicsmear -lStEpSimuJetEvent testClustering.cxx -o testClustering.exe

