#!/bin/bash

thiscut=$(echo $1 | sed "s?'??g")
scramdir=$2
outdir=$3
label=$(echo $thiscut | sed "s?[><&]?_?g")
thiscut=$(echo $thiscut | sed "s?&?\\\&?g")

echo "$thiscut"
echo "$scramdir" 
echo "$outdir"
echo "$label"

pwd
WD=$PWD

# cp -r $scramdir .
# ls
# cd CMSSW_7_4_7/src/MonoX
# scramv1 b ProjectRename
cd $scramdir/src/
eval `scramv1 runtime -sh`
# eval `/cvmfs/cms.cern.ch/common/scramv1 runtime -sh`
echo $CMSSW_BASE
which combine

cd $WD
cp -r $scramdir/src/MonoX .
cd MonoX
sed "s?XXXX?${thiscut}?g" configs/categories_config_vbf_tmpl.py > configs/categories_config_vbf_scan.py

python buildModel.py categories_config_vbf_scan
python runModel.py

cp combined_model.root datacards/ 
cd datacards
combine -M Asymptotic vbf_combined_noqcd.txt -n $label --run=blind --noFitAsimov 
#combine -M Asymptotic vbf_combined.txt -n $label --run=blind --noFitAsimov 
cp higgs*root $outdir 

cd $WD
rm -rf CMSSW_7_4_7
