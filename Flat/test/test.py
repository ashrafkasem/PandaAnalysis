#!/usr/bin/env python

from re import sub
from sys import argv,exit
from os import system,getenv
import json

torun = argv[1]
argv = []

import ROOT as root
from PandaCore.Tools.Load import *


Load('PandaAnalyzer')
print 'loaded'

skimmer = root.PandaAnalyzer()

print 'created'

#skimmer.firstEvent=0
#skimmer.lastEvent=10
skimmer.isData=False
skimmer.SetFlag('puppi',True)
skimmer.SetFlag('fatjet',True)
skimmer.SetFlag('firstGen',False)
skimmer.SetFlag('applyJSON',False)
if skimmer.isData and False:
    with open(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/certs/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt') as jsonFile:
        payload = json.load(jsonFile)
        for run,lumis in payload.iteritems():
            for l in lumis:
                skimmer.AddGoodLumiRange(int(run),l[0],l[1])
skimmer.processType = root.PandaAnalyzer.kTT
#        skimmer.SetPreselectionBit(root.PandaAnalyzer.kMonotop)
fin = root.TFile.Open(torun)

print torun
print fin

tree = fin.FindObjectAny("events")
hweights = fin.FindObjectAny("hSumW")
print tree,hweights

skimmer.SetDataDir(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/')
skimmer.SetOutputFile('testskim.root')
skimmer.Init(tree,hweights)

skimmer.Run()
print 'done running'
skimmer.Terminate()
print 'done terminating'

