#!/bin/bash 
#echo SingleTop Diboson WJets TTbar ZJets ZtoNuNu QCD MET SingleElectron SinglePhoton EWKZJets EWKZtoNuNu EWKWJets VBF_H125 GJets | xargs -n 1 -P 6 python merge.py
echo SingleTop Diboson WJets TTbar ZJets ZtoNuNu QCD MET SingleElectron | xargs -n 1 -P 6 python merge.py
