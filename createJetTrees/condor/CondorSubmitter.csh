#! /bin/env csh

set Exec = /direct/eic+u/bpage/epJets/createJetTrees/createJetTrees_heraJetShape_v3.exe
#set Exec = /direct/eic+u/bpage/epJets/createJetTrees/createJetTrees_algoRadTestPt_v3.exe
#set Exec = /direct/eic+u/bpage/epJets/createJetTrees/createJetTrees_aktR1.0Pt_v3.exe
#set Exec = /direct/eic+u/bpage/epJets/createJetTrees/createJetTrees_angRadPt_v3.exe

# make sure executable exists
#make $Exec || exit

####### Initialize condor file
echo ""  > CondorFile
echo "Universe     = vanilla" >> CondorFile
echo "Executable   = ${Exec}" >> CondorFile
echo "getenv = true" >> CondorFile
# echo "Notification = Complete" >> CondorFile
# echo "Notify_user  = kkauder@gmail.com"  >> CondorFile


# split into chunks
#set base = /gpfs/mnt/gpfs02/eic/bpage/extraSimu/PYTHIA/ep/TREES/*Q2=100.0-500.0*
#set base = /gpfs/mnt/gpfs02/eic/bpage/extraSimu/PYTHIA/ep/TREES/*5Mevents*
#set base = /gpfs/mnt/gpfs02/eic/bpage/extraSimu/PYTHIA/ep/TREES/*10x50*Q2=1.0-10.0*
#set base = /eicdata/eic0009/PYTHIA/ep/TREES/RealPhotons/Q2_10-5_1/pythia.ep.20x250.5Mevents*
set base = /gpfs/mnt/gpfs02/eic/bpage/extraSimu/PYTHIA/ep/TREES/HERA/*
#set base = /gpfs/mnt/gpfs02/eic/bpage/extraSimu/PYTHIA/ep/TREES/ISRFSR/*.root
#set base = /gpfs/mnt/gpfs02/eic/bpage/extraSimu/PYTHIA/ep/TREES/HADISRFSR/*NOISRFSR*.root
#set base = /gpfs/mnt/gpfs02/eic/bpage/extraSimu/PYTHIA/ep/TREES/*20x250.1Mevents*Q2=0.00001-1.0*

# Output Directory
#set Output = /direct/eic+u/bpage/epJets/jetTrees_v3/kinematics/condor
#set Output = /gpfs/mnt/gpfs02/eic/bpage/test
#set Output = /gpfs/mnt/gpfs02/eic/bpage/algoRadiusTest/Q2_1-10/extra/10x50
#set Output = /gpfs/mnt/gpfs02/eic/bpage/algoRadiusTest/Q2_10-100/extra
#set Output = /gpfs/mnt/gpfs02/eic/bpage/algoRadiusTest/Q2_100-500/extra
#set Output = /gpfs/mnt/gpfs02/eic/bpage/algoRadiusTest/Q2_10-5_1
#set Output = /gpfs/mnt/gpfs02/eic/bpage/algoRadiusTest/Root6
set Output = /gpfs/mnt/gpfs02/eic/bpage/heraJetShape
#set Output = /gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10Pt_v3/Root6
#set Output = /gpfs/mnt/gpfs02/eic/bpage/jetTrees_angRadPt

foreach input ( ${base}* )
    # arguments
    set OutBasePre=`basename $input | sed 's/.root//g'`
    set OutBase=`echo $OutBasePre | sed 's/pythia//'`
    set OutNameA    = ${Output}/jetTree${OutBase}.A.root
    set OutNameB    = ${Output}/jetTree${OutBase}.B.root
    set OutNameC    = ${Output}/jetTree${OutBase}.C.root
    set OutNameD    = ${Output}/jetTree${OutBase}.D.root
    set OutName    = ${Output}/jetTree${OutBase}.T.root
    set Files      = ${input}
    
    # Logfiles.
    set LogFileA    = ${Output}/logs/jetTree${OutBase}.A.out
    set LogFileB    = ${Output}/logs/jetTree${OutBase}.B.out
    set LogFileC    = ${Output}/logs/jetTree${OutBase}.C.out
    set LogFileD    = ${Output}/logs/jetTree${OutBase}.D.out
    set LogFile    = ${Output}/logs/jetTree${OutBase}.T.out
    set ErrFileA    = ${Output}/logs/jetTree${OutBase}.A.err
    set ErrFileB    = ${Output}/logs/jetTree${OutBase}.B.err
    set ErrFileC    = ${Output}/logs/jetTree${OutBase}.C.err
    set ErrFileD    = ${Output}/logs/jetTree${OutBase}.D.err
    set ErrFile    = ${Output}/logs/jetTree${OutBase}.T.err

    ### hand to condor
    #set Args = ( -o $OutName -i $Files -breit $breit )
    #set ArgsA = ( 250000 0 250000 $Files $OutNameA )
    set ArgsA = ( 1000000 0 250000 $Files $OutNameA )
    set ArgsB = ( 1000000 250000 500000 $Files $OutNameB )
    set ArgsC = ( 1000000 500000 750000 $Files $OutNameC )
    set ArgsD = ( 1000000 750000 1000000 $Files $OutNameD )
    set Args = ( 5000000 4750000 5000000 $Files $OutName )
    echo "" >> CondorFile
    echo "Output       = ${LogFileD}" >> CondorFile
    echo "Error        = ${ErrFileD}" >> CondorFile
    echo "Arguments    = ${ArgsD}" >> CondorFile
    echo "Queue" >> CondorFile   

    echo Submitting:
    echo $Exec $ArgsD
    echo "Logging output to " $LogFileD
    echo "Logging errors to " $ErrFileD
    echo
end
#condor_submit CondorFile


