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
   
    skimmer.firstEvent=0
    skimmer.lastEvent=100
    skimmer.isData=False
    skimmer.doHeatMap=False
    skimmer.doECF=True
    skimmer.doQjets=False
    skimmer.doKinFit=False
    if 'TT' in fullPath:
      skimmer.processType = root.Analyzer.kTop
    elif 'WW' in fullPath:
      skimmer.processType = root.Analyzer.kV
    else:
      skimmer.processType = root.Analyzer.kQCD
    skimmer.processType = root.Analyzer.kQCD
    eosPath = 'root://eoscms//eos/cms/store/user/%s'%(getenv('USER'))
    cernboxPath = 'root://eosuser//eos/user/%s/%s'%(getenv('USER')[0],getenv('USER'))
    cernboxBPath = 'root://eosuser//eos/user/b/bmaier'

    fullPath = sub(r'\${CERNBOXB}',cernboxBPath,sub(r'\${CERNBOX}',cernboxPath,sub(r'\${EOS}',eosPath,fullPath)))
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

  fn(argv[1]) 
