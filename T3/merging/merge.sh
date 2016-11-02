#!/bin/bash 
#echo MET Diboson WJets TTbar SingleTop QCD SinglePhoton GJets  | xargs -n 1 -P 12 python merge.py
echo MET SinglePhoton GJets WJets TTbar SingleTop QCD Diboson ZtoNuNu SingleElectron | xargs -n 1 -P 8 python merge.py
