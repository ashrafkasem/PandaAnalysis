#!/bin/bash

fittingdir=$1
scramdir=$2
model=$3
mV=$4
mChi=$5

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
which combine

cd $WD
cp -r $scramdir/src/MonoXFit_CSV .
cd MonoXFit_CSV/datacards

python scanbatch.py tmpl.txt --mMed $mV --mChi $mChi --infile $fittingdir/fittingForest.root $model

#cp signalmodel.root $fittingdir/scans/signal_${mV}_${mChi}.root
cp higgs*root $fittingdir/scans
#cp scan_*txt $fittingdir/scans

cd $WD
#cp -r MonoXFit_CSV $fittingdir/scans/fit_${mV}_${mChi}_${model}
rm -rf MonoXFit_CSV
