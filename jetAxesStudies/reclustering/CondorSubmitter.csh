#! /bin/env csh

#set Exec = ./bin/RunTest
set Exec = /direct/eic+u/bpage/epJets/jetTrees_v3/jetAxesStudy/analyzeReclusterWTA_v0.exe

# make sure executable exists
#make $Exec || exit

####### Initialize condor file
echo ""  > CondorFile
echo "Universe     = vanilla" >> CondorFile
echo "Executable   = ${Exec}" >> CondorFile
echo "getenv = true" >> CondorFile
# echo "Notification = Complete" >> CondorFile
# echo "Notify_user  = kkauder@gmail.com"  >> CondorFile

#set breit = 1

#switch ($breit)
#    case 0 : 
#    set NameBase=Lab
#    breaksw
#    case 1 : 
#    set NameBase=Breit
#    breaksw
#    default : 
#    echo "unknown breit setting $breit"
#    exit -1
#endsw


# split into chunks
#set base = Data/extraSimu/PYTHIA/ep/TREES/pythia.ep.20x250.*Q2=10.0-100*root
#set base = Data/extraSimu/PYTHIA/ep/TREES/pythia.ep.20x250.*root
#set base = /gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10_v3/Q2_10-100/particle/jetTree*kT=1.0_*.root
set base = /gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10Pt_v3/Root6/Q2_10-100/particle/jetTree.ep.20x250.*.root
#set base = /gpfs/mnt/gpfs02/eic/bpage/jetTrees_akt10Pt_v3/Root6/*.root

# Output Directory
#set Output = /direct/eic+u/bpage/epJets/jetTrees_v3/kinematics/condor
set Output = /gpfs/mnt/gpfs02/eic/bpage/analysis

foreach input ( ${base}* )
    # arguments
    set OutBase=`basename $input | sed 's/.root//g'`
    set OutName    = ${Output}/out_${OutBase}.hist.root
    set Files      = ${input}
    
    # Logfiles.
    set LogFile    = ${Output}/logs/out_${OutBase}.out
    set ErrFile    = ${Output}/logs/out_${OutBase}.err

    ### hand to condor
    #set Args = ( -o $OutName -i $Files -breit $breit )
    set Args = ( 0 $Files $OutName )
    echo "" >> CondorFile
    echo "Output       = ${LogFile}" >> CondorFile
    echo "Error        = ${ErrFile}" >> CondorFile
    echo "Arguments    = ${Args}" >> CondorFile
    echo "Queue" >> CondorFile   

    echo Submitting:
    echo $Exec $Args
    echo "Logging output to " $LogFile
    echo "Logging errors to " $ErrFile
    echo
end
#condor_submit CondorFile


