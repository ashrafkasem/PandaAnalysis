#!/bin/bash 
#echo MET Diboson WJets TTbar SingleTop QCD SinglePhoton GJets  | xargs -n 1 -P 12 python merge.py
echo MET SinglePhoton SingleElectron TTbar ZtoNuNu ZJets GJets WJets SingleTop QCD Diboson TTbar_isrup TTbar_isrdown TTbar_tuneup TTbar_tunedown | xargs -n 1 -P 10 python merge.py
#echo MET SinglePhoton TTbar ZtoNuNu ZJets GJets WJets SingleTop QCD Diboson TTbar_isrup TTbar_isrdown TTbar_FXFX TTbar_tuneup TTbar_tunedown | xargs -n 1 -P 10 python merge.py
#echo MET SingleElectron SinglePhoton | xargs -n 1 -P 5 python merge.py
