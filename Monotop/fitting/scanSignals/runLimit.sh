#!/bin/bash

fittingdir=$1
scramdir=$2
model=$3
mParams=$4
couplings=$5

echo "ARGS $@"

echo -n "PWD "
pwd
WD=$PWD

#cp -r $scramdir .
#ls
#cd CMSSW_7_4_7/src/MonoX
#scramv1 b ProjectRename
cd $scramdir/src/
eval `scramv1 runtime -sh`
# eval `/cvmfs/cms.cern.ch/common/scramv1 runtime -sh`
echo $CMSSW_BASE
echo -n "COMBINE "
which combine

cd $WD
cp -r $scramdir/src/MonoXFit_CSV .
cd MonoXFit_CSV/datacards

if [[ "${couplings}" == "" ]]; then 
    python scanbatch.py correlated_tmpl.txt --mParams $mParams --infile $fittingdir/fittingForest.root --model $model
else
    python scanbatch.py correlated_tmpl.txt --mParams $mParams --infile $fittingdir/signals/fittingForest_signal_vector_${couplings}_nlo.root --model $model
fi

#cp signalmodel.root $fittingdir/scans/signal_${mV}_${mChi}.root
if [[ "${couplings}" == "" ]]; then 
    mkdir -p $fittingdir/scans/$model/nominal/
    cp higgs*root $fittingdir/scans/$model/nominal/
else
    mkdir -p $fittingdir/scans/$model/${couplings}/
    cp higgs*root $fittingdir/scans/$model/${couplings}/
fi
#cp scan_*txt $fittingdir/scans

cd $WD
#cp -r MonoXFit_CSV $fittingdir/scans/fit_${mV}_${mChi}_${model}
rm -rf MonoXFit_CSV
