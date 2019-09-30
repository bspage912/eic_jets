rootcint -f out.cxx -c -p StEpSimuJetEvent.h StEpSimuJetDef.h StEpSimuJet.h StEpSimuJetParticle.h StEpSimuParticle.h StEpSimuSubJet.h LinkDef.h

g++ -Wall -ggdb -fPIC `root-config --cflags --libs` -c StEpSimuJetEvent.cxx StEpSimuJetDef.cxx StEpSimuJet.cxx StEpSimuJetParticle.cxx StEpSimuParticle.cxx StEpSimuSubJet.cxx out.cxx

c++ -Wall -ggdb -shared `root-config --cflags --libs` -o libStEpSimuJetEvent.so StEpSimuJetEvent.o StEpSimuJetDef.o StEpSimuJet.o StEpSimuJetParticle.o StEpSimuParticle.o StEpSimuSubJet.o out.o

rm *.o
rm out.*

#c++ -Wall -ggdb `root-config --cflags --libs` -I/eic/u/bpage/epJets/test1 -L/eic/u/bpage/epJets/test1 -lrootWrapper -lfastjetWrapper runFastJetWrapperRoot.cxx -o runRoot.exe