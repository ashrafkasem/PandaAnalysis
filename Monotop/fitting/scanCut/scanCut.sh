#!/bin/bash

thiscut=$(echo $1 | sed "s?'??g")
scramdir=$2
outdir=$3
label="_$(echo $thiscut | sed "s?[\(\)_]??g" | sed "s?[><&]?_?g" )_"
#thiscut=$(echo $thiscut | sed "s?&?\\\&?g")

echo "$thiscut"
echo "$scramdir" 
echo "$outdir"
echo "$label"

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
cd MonoXFit_CSV
#cd $CMSSW_BASE/src/MonoXFit_CSV/
#sed "s?XXXX?${thiscut}?g" configs/categories_tmpl.py > configs/categories_scan.py

#python buildModel.py test

cp ${outdir}/hists/mono-x_${thiscut}.root mono-x.root
python runModel.py

cd datacards
combine -M Asymptotic combined_wtop.txt -n $label --run=blind --noFitAsimov 
cp higgs*root $outdir 

#cp -r ../MonoXFit_CSV/ $outdir/fit_$label

cd $WD
