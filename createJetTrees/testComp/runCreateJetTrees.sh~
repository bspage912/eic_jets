#!/bin/sh

number=51
limit=100

# OLD LOCAL
#path=/eicdata/eic0009/PYTHIA/ep/TREES/Q2-KT-BINS
#dest=/eicdata/eic0009/bpage/jetTrees_akt10_v3/Q2_10-100
#dest=/eicdata/eic0009/bpage/jetTrees_akt10_v3/Q2_1-10
#dest=/eicdata/eic0009/bpage/jetTrees_akt10_v3/Q2_100-500
#dest=/eicdata/eic0009/bpage/eektTest/Q2_10-100


#path=/eicdata/eic0009/PYTHIA/ep/TREES/RealPhotons
#dest=/eicdata/eic0009/bpage/jetTrees_akt10_v3/lowQ2
#dest=/eicdata/eic0009/bpage/jetTrees_akt10_v2/hera
#dest=/eicdata/eic0009/bpage/eektTest/lowQ2


#path=/eicdata/eic0009/PYTHIA/ep/TREES/
#dest=/eicdata/eic0009/bpage/jetTrees_akt10_v3/

#path=/eicdata/eic0009/PYTHIA/ep/TREES/Q2-KT-BINS
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_10-100/particle_noNeutral

# EXTRA SIMU
path=/gpfs/mnt/gpfs02/eic/bpage/extraSimu/PYTHIA/ep/TREES
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_100-500/particle/extra
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_10-100/particle/extra
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_1-10/particle/extra
dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_01-1/particle/extra

# SMEARING
#path=/eicdata/eic0009/PYTHIA/ep/TREES/Q2-KT-BINS
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_10-100/smeared_beast
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_10-100/trackEff95
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_10-100/trackEff90
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_10-100/keepAllNeutrals
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_10-100/keepAllNeutrals_zeus
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_10-100/smeared_beast_noNeutral
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_1-10/smeared_beast
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_1-10/trackEff95
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_1-10/trackEff90
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_1-10/keepAllNeutrals
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3
#dest=/direct/eic+u/bpage/epJets/createJetTrees

# Algo Radius Test
#path=/eicdata/eic0009/PYTHIA/ep/TREES/Q2-KT-BINS
#dest=/gpfs/mnt/gpfs02/eic/bpage/algoRadiusTest/Q2_10-100
#dest=/gpfs/mnt/gpfs02/eic/bpage/algoRadiusTest/Q2_1-10
#dest=/gpfs/mnt/gpfs02/eic/bpage/algoRadiusTest/Q2_100-500
#dest=/gpfs/mnt/gpfs02/eic/bpage/algoRadiusTest/Q2_500-1000

# Parton Jets
#path=/eicdata/eic0009/PYTHIA/ep/TREES/RealPhotons/Q2_10-5_1
#dest=/gpfs/mnt/gpfs02/eic/bpage/jetTrees_partonRad


# Make Temporary Local Directory
mkdir /data1/pagetmp


while [ "$number" -le "$limit" ]
do 

    echo $number

    date

    # Low Q2 Trees
    #./createJetTrees_aktR1.0_v3.exe 5000000 0 250000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 250000 500000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 500000 750000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 750000 1000000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.D.root

    #./createJetTrees_aktR1.0_v3.exe 5000000 1000000 1250000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.E.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 1250000 1500000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.F.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 1500000 1750000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.G.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 1750000 2000000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.H.root

    #./createJetTrees_aktR1.0_v3.exe 5000000 2000000 2250000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.I.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 2250000 2500000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.J.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 2500000 2750000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.K.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 2750000 3000000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.L.root

    #./createJetTrees_aktR1.0_v3.exe 5000000 3000000 3250000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.M.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 3250000 3500000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.N.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 3500000 3750000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.O.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 3750000 4000000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.P.root

    #./createJetTrees_aktR1.0_v3.exe 5000000 4000000 4250000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.Q.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 4250000 4500000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.R.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 4500000 4750000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.S.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 4750000 5000000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.T.root

    # Q2 = 0.01 - 0.1 kT = 1.0
    #./createJetTrees_aktR1.0_v3.exe 1000000 0 250000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=0.01-0.1.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.01-0.1.kT=1.0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 250000 500000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=0.01-0.1.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.01-0.1.kT=1.0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 500000 750000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=0.01-0.1.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.01-0.1.kT=1.0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 750000 1000000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=0.01-0.1.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.01-0.1.kT=1.0_$number.D.root

    # Q2 = 0.1 - 1.0 kT = 1.0
    ./createJetTrees_aktR1.0_v3.exe 1000000 0 250000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.A.root
    ./createJetTrees_aktR1.0_v3.exe 1000000 250000 500000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.B.root
    ./createJetTrees_aktR1.0_v3.exe 1000000 500000 750000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.C.root
    ./createJetTrees_aktR1.0_v3.exe 1000000 750000 1000000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.D.root

    #./createJetTrees_aktR1.0_v3.exe 1000000 0 250000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 250000 500000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 500000 750000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 750000 1000000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=0.1-1.0.kT=1.0_$number.D.root

    # Q2 = 1 - 10 kT = 2.0
    #./createJetTrees_aktR1.0_v3.exe 1000000 0 250000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=2.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=2.0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 250000 500000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=2.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=2.0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 500000 750000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=2.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=2.0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 750000 1000000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=2.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=2.0_$number.D.root

    # Q2 = 1 - 10 kT = 1.0
    #./createJetTrees_aktR1.0_v3.exe 1000000 0 250000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 250000 500000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 500000 750000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 750000 1000000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.D.root

    #./createJetTreesSmeared_aktR1.0_v3.exe 1000000 0 250000 1 1.00 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.smear.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.A.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 1000000 250000 500000 1 1.00 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.smear.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.B.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 1000000 500000 750000 1 1.00 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.smear.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.C.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 1000000 750000 1000000 1 1.00 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.smear.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.D.smear.root

    #./createJetTrees_algoRadTest_v3.exe 1000000 0 250000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.A.root
    #./createJetTrees_algoRadTest_v3.exe 1000000 250000 500000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.B.root
    #./createJetTrees_algoRadTest_v3.exe 1000000 500000 750000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.C.root
    #./createJetTrees_algoRadTest_v3.exe 1000000 750000 1000000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.D.root

    #./createJetTrees_aktR1.0_v3.exe 1000000 0 250000 $path/pythia.ep.10x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.10x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 250000 500000 $path/pythia.ep.10x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.10x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 500000 750000 $path/pythia.ep.10x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.10x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 750000 1000000 $path/pythia.ep.10x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.10x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.D.root

    #./createJetTrees_aktR1.0_v3.exe 1000000 0 250000 $path/pythia.ep.20x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.20x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 250000 500000 $path/pythia.ep.20x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.20x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 500000 750000 $path/pythia.ep.20x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.20x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 750000 1000000 $path/pythia.ep.20x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.20x100.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.D.root

    #./createJetTrees_aktR1.0_v3.exe 1000000 0 250000 $path/pythia.ep.10x40.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.10x40.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 250000 500000 $path/pythia.ep.10x40.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.10x40.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 500000 750000 $path/pythia.ep.10x40.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.10x40.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 750000 1000000 $path/pythia.ep.10x40.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.root jetTree.ep.10x40.1Mevents.RadCor=0.Q2=1.0-10.0.kT=1.0_$number.D.root

    # Q2 = 10 - 100 kT = 2.0
    #./createJetTrees_aktR1.0_v3.exe 1000000 0 250000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=2.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=2.0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 250000 500000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=2.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=2.0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 500000 750000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=2.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=2.0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 750000 1000000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=2.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=2.0_$number.D.root

    # Q2 = 10 - 100 kT = 1.0
    #./createJetTrees_aktR1.0_v3.exe 1000000 0 250000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 250000 500000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 500000 750000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 750000 1000000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.D.root

    #./createJetTreesSmeared_aktR1.0_v3.exe 1000000 0 250000 1 1.00 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.A.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 1000000 250000 500000 1 1.00 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.B.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 1000000 500000 750000 1 1.00 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.C.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 1000000 750000 1000000 1 1.00 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.D.smear.root

    #./createJetTrees_algoRadTest_v3.exe 1000000 0 250000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.A.root
    #./createJetTrees_algoRadTest_v3.exe 1000000 250000 500000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.B.root
    #./createJetTrees_algoRadTest_v3.exe 1000000 500000 750000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.C.root
    #./createJetTrees_algoRadTest_v3.exe 1000000 750000 1000000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.D.root

    #./createJetTrees_aktR1.0_v3.exe 1000000 0 250000 $path/pythia.ep.10x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.10x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 250000 500000 $path/pythia.ep.10x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.10x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 500000 750000 $path/pythia.ep.10x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.10x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 750000 1000000 $path/pythia.ep.10x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.10x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.D.root

    #./createJetTrees_aktR1.0_v3.exe 1000000 0 250000 $path/pythia.ep.20x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 250000 500000 $path/pythia.ep.20x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 500000 750000 $path/pythia.ep.20x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 750000 1000000 $path/pythia.ep.20x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x100.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.D.root

    #./createJetTrees_aktR1.0_v3.exe 1000000 0 250000 $path/pythia.ep.10x40.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.10x40.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 250000 500000 $path/pythia.ep.10x40.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.10x40.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 500000 750000 $path/pythia.ep.10x40.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.10x40.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 1000000 750000 1000000 $path/pythia.ep.10x40.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.10x40.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.D.root

    # Q2 = 10 - 100 kT = 1.0 5M Event Files
    #./createJetTrees_aktR1.0_v3.exe 5000000 0 250000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 250000 500000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 500000 750000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 750000 1000000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.D.root

    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 0 250000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.A.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 250000 500000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.B.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 500000 750000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.C.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 750000 1000000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.D.smear.root

    #./createJetTrees_aktR1.0_v3.exe 5000000 1000000 1250000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.E.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 1250000 1500000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.F.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 1500000 1750000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.G.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 1750000 2000000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.H.root

    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 1000000 1250000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.E.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 1250000 1500000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.F.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 1500000 1750000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.G.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 1750000 2000000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.H.smear.root

    #./createJetTrees_aktR1.0_v3.exe 5000000 2000000 2250000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.I.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 2250000 2500000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.J.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 2500000 2750000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.K.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 2750000 3000000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.L.root

    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 2000000 2250000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.I.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 2250000 2500000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.J.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 2500000 2750000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.K.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 2750000 3000000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.L.smear.root

    #./createJetTrees_aktR1.0_v3.exe 5000000 3000000 3250000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.M.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 3250000 3500000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.N.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 3500000 3750000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.O.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 3750000 4000000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.P.root

    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 3000000 3250000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.M.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 3250000 3500000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.N.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 3500000 3750000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.O.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 3750000 4000000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.P.smear.root

    #./createJetTrees_aktR1.0_v3.exe 5000000 4000000 4250000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.Q.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 4250000 4500000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.R.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 4500000 4750000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.S.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 4750000 5000000 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.T.root

    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 4000000 4250000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.Q.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 4250000 4500000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.R.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 4500000 4750000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.S.smear.root
    #./createJetTreesSmeared_aktR1.0_v3.exe 5000000 4750000 5000000 1 1.00 $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root $path/pythia.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.smearB.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.T.smear.root

    # Q2 = 100 - 500 kT = 1.0
    #./createJetTrees_aktR1.0_v3.exe 250000 0 250000 $path/pythia.ep.20x250.250Kevents.RadCor=0.Q2=100.0-500.0.kT=1.0_$number.root /data1/pagetmp/jetTree.ep.20x250.250Kevents.RadCor=0.Q2=100.0-500.0.kT=1.0_$number.A.root

    #./createJetTrees_aktR1.0_v3.exe 250000 0 250000 $path/pythia.ep.20x250.250Kevents.RadCor=0.Q2=100.0-500.0.kT=1.0_$number.root jetTree.ep.20x250.250Kevents.RadCor=0.Q2=100.0-500.0.kT=1.0_$number.A.root

    #./createJetTrees_algoRadTest_v3.exe 250000 0 250000 $path/pythia.ep.20x250.250Kevents.RadCor=0.Q2=100.0-500.0.kT=1.0_$number.root jetTree.ep.20x250.250Kevents.RadCor=0.Q2=100.0-500.0.kT=1.0_$number.A.root

    # Q2 = 500 - 1000 kT = 1.0
    #./createJetTrees_algoRadTest_v3.exe 50000 0 50000 $path/pythia.ep.20x250.50Kevents.RadCor=0.Q2=500.0-1000.0.kT=1.0_$number.root jetTree.ep.20x250.50Kevents.RadCor=0.Q2=500.0-1000.0.kT=1.0_$number.A.root

    # No Q2 Cut
    #./createJetTrees_aktR1.0_v3.exe 5000000 0 250000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.root jetTree.ep.20x250.5Mevents.RadCor=0_$number.A.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 250000 500000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.root jetTree.ep.20x250.5Mevents.RadCor=0_$number.B.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 500000 750000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.root jetTree.ep.20x250.5Mevents.RadCor=0_$number.C.root
    #./createJetTrees_aktR1.0_v3.exe 5000000 750000 1000000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.root jetTree.ep.20x250.5Mevents.RadCor=0_$number.D.root

    # EEKT Test Q2 = 10 - 100
    #./createJetTrees_eektTest.exe 1000000 0 250000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.A.root
    #./createJetTrees_eektTest.exe 1000000 250000 500000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.B.root
    #./createJetTrees_eektTest.exe 1000000 500000 750000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.C.root
    #./createJetTrees_eektTest.exe 1000000 750000 1000000 $path/pythia.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.root jetTree.ep.20x250.1Mevents.RadCor=0.Q2=10.0-100.0.kT=1.0_$number.D.root

    # EEKT Test Low Q2 Trees
    #./createJetTrees_eektTest.exe 5000000 0 250000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.A.root
    #./createJetTrees_eektTest.exe 5000000 250000 500000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.B.root
    #./createJetTrees_eektTest.exe 5000000 500000 750000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.C.root
    #./createJetTrees_eektTest.exe 5000000 750000 1000000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.D.root

    #./createJetTrees_eektTest.exe 5000000 1000000 1250000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.E.root
    #./createJetTrees_eektTest.exe 5000000 1250000 1500000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.F.root
    #./createJetTrees_eektTest.exe 5000000 1500000 1750000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.G.root
    #./createJetTrees_eektTest.exe 5000000 1750000 2000000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.H.root

    #./createJetTrees_eektTest.exe 5000000 2000000 2250000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.I.root
    #./createJetTrees_eektTest.exe 5000000 2250000 2500000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.J.root
    #./createJetTrees_eektTest.exe 5000000 2500000 2750000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.K.root
    #./createJetTrees_eektTest.exe 5000000 2750000 3000000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.L.root

    #./createJetTrees_eektTest.exe 5000000 3000000 3250000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.M.root
    #./createJetTrees_eektTest.exe 5000000 3250000 3500000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.N.root
    #./createJetTrees_eektTest.exe 5000000 3500000 3750000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.O.root
    #./createJetTrees_eektTest.exe 5000000 3750000 4000000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.P.root

    #./createJetTrees_eektTest.exe 5000000 4000000 4250000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.Q.root
    #./createJetTrees_eektTest.exe 5000000 4250000 4500000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.R.root
    #./createJetTrees_eektTest.exe 5000000 4500000 4750000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.S.root
    #./createJetTrees_eektTest.exe 5000000 4750000 5000000 $path/pythia.ep.20x250.5Mevents.1.RadCor=0.Q2-0.-10_SAS_$number.root jetTree.ep.20x250.5Mevents.RadCor=0.Q2=0.0-10.0_$number.T.root

    # Hera Trees
    #./createJetTrees_aktR1.0_v2.exe 5000000 0 250000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.A.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 250000 500000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.B.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 500000 750000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.C.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 750000 1000000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.D.root

    #./createJetTrees_aktR1.0_v2.exe 5000000 1000000 1250000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.E.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 1250000 1500000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.F.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 1500000 1750000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.G.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 1750000 2000000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.H.root

    #./createJetTrees_aktR1.0_v2.exe 5000000 2000000 2250000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.I.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 2250000 2500000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.J.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 2500000 2750000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.K.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 2750000 3000000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.L.root

    #./createJetTrees_aktR1.0_v2.exe 5000000 3000000 3250000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.M.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 3250000 3500000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.N.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 3500000 3750000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.O.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 3750000 4000000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.P.root

    #./createJetTrees_aktR1.0_v2.exe 5000000 4000000 4250000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.Q.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 4250000 4500000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.R.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 4500000 4750000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.S.root
    #./createJetTrees_aktR1.0_v2.exe 5000000 4750000 5000000 $path/HERA1.root jetTree.ep.27x820.5Mevents_$number.T.root


    #./createJetTrees_partonRad_v3.exe 5000000 0 250000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.A.root
    #./createJetTrees_partonRad_v3.exe 5000000 250000 500000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.B.root
    #./createJetTrees_partonRad_v3.exe 5000000 500000 750000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.C.root
    #./createJetTrees_partonRad_v3.exe 5000000 750000 1000000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.D.root

    #./createJetTrees_partonRad_v3.exe 5000000 1000000 1250000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.E.root
    #./createJetTrees_partonRad_v3.exe 5000000 1250000 1500000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.F.root
    #./createJetTrees_partonRad_v3.exe 5000000 1500000 1750000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.G.root
    #./createJetTrees_partonRad_v3.exe 5000000 1750000 2000000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.H.root

    #./createJetTrees_partonRad_v3.exe 5000000 2000000 2250000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.I.root
    #./createJetTrees_partonRad_v3.exe 5000000 2250000 2500000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.J.root
    #./createJetTrees_partonRad_v3.exe 5000000 2500000 2750000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.K.root
    #./createJetTrees_partonRad_v3.exe 5000000 2750000 3000000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.L.root

    #./createJetTrees_partonRad_v3.exe 5000000 3000000 3250000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.M.root
    #./createJetTrees_partonRad_v3.exe 5000000 3250000 3500000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.N.root
    #./createJetTrees_partonRad_v3.exe 5000000 3500000 3750000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.O.root
    #./createJetTrees_partonRad_v3.exe 5000000 3750000 4000000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.P.root

    #./createJetTrees_partonRad_v3.exe 5000000 4000000 4250000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.Q.root
    #./createJetTrees_partonRad_v3.exe 5000000 4250000 4500000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.R.root
    #./createJetTrees_partonRad_v3.exe 5000000 4500000 4750000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.S.root
    #./createJetTrees_partonRad_v3.exe 5000000 4750000 5000000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root /data1/pagetmp/jetTree.ep.20x250.5Mevents.RadCor=0_$number.T.root

    #./createJetTrees_partonRad_v3.exe 250000 0 250000 $path/pythia.ep.20x250.5Mevents.$number.RadCor=0.Q2-0.-1_GRV.root jetTree.ep.20x250.5Mevents.RadCor=0_$number.A.root

    
    date

    mv /data1/pagetmp/*.root $dest

    date

    let number++

done

# Remove Temp Directory
rmdir /data1/pagetmp


