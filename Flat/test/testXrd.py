#!/usr/bin/env python

from ROOT import gSystem,gROOT
import ROOT as root
from re import sub
from sys import argv,exit
from os import system,getenv
from PandaCore.Tools.Load import *

if __name__ == "__main__":

  Load('SCRAMJetAnalyzer','Analyzer')
#  gROOT.LoadMacro("${CMSSW_BASE}/src/SCRAMJet/Analyzer/interface/Analyzer.h")
#  gSystem.Load('libSCRAMJetAnalyzer.so')
#  gSystem.Load('libSCRAMJetObjects.so')
  
  def fn(fullPath):

    skimmer = root.Analyzer()
   
    skimmer.maxEvents = 10
    skimmer.isData=False
    skimmer.doHeatMap=True
    skimmer.doECF=True
    skimmer.doQjets=True
    skimmer.doKinFit=False
    if 'TT' in fullPath:
      skimmer.processType = root.Analyzer.kTop
    elif 'WW' in fullPath:
      skimmer.processType = root.Analyzer.kV
    else:
      skimmer.processType = root.Analyzer.kQCD
    fin = root.TFile.Open(fullPath)

    print fullPath
    print fin

    tree = fin.FindObjectAny("events")
    print tree

    skimmer.SetOutputFile('/tmp/%s/testskim.root'%getenv('USER'))
    skimmer.Init(tree)
    skimmer.AddFatJetFromTree("puppiCA15","puppiCA15",root.Analyzer.kPuppi,1.5,root.Analyzer.kCA)

    skimmer.Run()
    print 'done running'
    skimmer.Terminate()
    print 'done terminating'

  fn('root://eoscms.cern.ch//store/user/snarayan/scramjet/ZprimeToTTJet_M-1000_TuneCUETP8M1_13TeV-amcatnlo-pythia8/ZprimeToTTJet_M-1000_TuneCUETP8M1_13TeV-amcatnlo-pythia8/160817_205018/0000/scramjet_1.root') 
