# g++ -Wall -ggdb `root-config --cflags --libs` \
#     `fastjet-config --cxxflags --libs --plugins` \
#     -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -lStEpSimuJetEvent \
#     -lRecursiveTools \
#     analyzeReclusterWTA_v0.cxx -o analyzeReclusterWTA_v0.exe

g++  -o analyzeReclusterWTA_v0.exe analyzeReclusterWTA_v0.cxx \
    -Wall -ggdb `root-config --cflags --libs` \
    -I/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -L/eic/u/bpage/epJets/createJetTrees/StEpSimuJetMaker_Root6 -lStEpSimuJetEvent \
    `fastjet-config --cxxflags --libs --plugins` \
    -lRecursiveTools 

#g++ -Wall -ggdb `root-config --cflags --libs` -I/direct/eic+data/foloic/wrk/150818/forXiaoxuan/createJetTrees/StEpSimuJetMaker -L/direct/eic+data/foloic/wrk/150818/forXiaoxuan/createJetTrees/StEpSimuJetMaker -lStEpSimuJetEvent analyzeJetKinematics_xiaoxuan_v1.cxx -o analyzeJetKinematics_xiaoxuan_v1.exe
#    -std=c++11 






