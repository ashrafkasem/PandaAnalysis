#!/bin/bash 
#echo MET Diboson WJets TTbar SingleTop QCD SinglePhoton GJets  | xargs -n 1 -P 12 python merge.py
echo MET SinglePhoton GJets WJets TTbar SingleTop QCD Diboson | xargs -n 1 -P 6 python merge.py
