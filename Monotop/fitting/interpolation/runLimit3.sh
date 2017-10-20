#!/bin/bash

scramdir=$1

echo "ARGS $@"

echo -n "PWD "
pwd
WD=$PWD

#cp -r $scramdir .
#ls
#cd CMSSW_7_4_7/src/MonoX
#scramv1 b ProjectRename
cd $scramdir/src/
eval `/cvmfs/cms.cern.ch/common/scramv1 runtime -sh`
echo $CMSSW_BASE
echo -n "COMBINE "
which combine

cd $WD
cp -r $scramdir/src/MonoXFit_CSV .
cd MonoXFit_CSV/datacards

python scanbatch3.py "${@:2}"
ret=$?

cd $WD
#cp -r MonoXFit_CSV $fittingdir/scans/fit_${mV}_${mChi}_${model}
rm -rf MonoXFit_CSV

exit $ret
