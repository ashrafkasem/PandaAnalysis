#!/bin/bash 
echo Diboson WJets TTbar ZJets ZtoNuNu QCD MET SingleElectron SinglePhoton EWKZJets EWKZtoNuNu EWKWJets VBF_H125 GJets | xargs -n 1 -P 12 python merge.py
