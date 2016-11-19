#!/bin/bash 
#echo MET Diboson WJets TTbar SingleTop QCD SinglePhoton GJets  | xargs -n 1 -P 12 python merge.py
echo MET SinglePhoton GJets WJets TTbar TTbar_Herwig SingleTop QCD Diboson ZtoNuNu ZJets SingleElectron monotop-nr-v3-1700-100_med-1700_dm-100 | xargs -n 1 -P 20 python merge.py
