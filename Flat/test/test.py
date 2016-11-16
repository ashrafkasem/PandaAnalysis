#!/usr/bin/env python

from ROOT import gSystem,gROOT
import ROOT as root
from re import sub
from sys import argv,exit
from os import system,getenv
from PandaCore.Tools.Load import *

if __name__ == "__main__":

  Load('PandaAnalysisFlat','PandaAnalyzer')
  
  def fn(fullPath):

    skimmer = root.PandaAnalyzer()
   
    skimmer.firstEvent=0
    skimmer.lastEvent=100
    skimmer.isData=False
    skimmer.applyJson=False
    skimmer.doFatjet=False
    skimmer.SetPreselectionBit(root.PandaAnalyzer.kMonojet)
#    skimmer.SetPreselectionBit(root.PandaAnalyzer.kMonotop)
#    skimmer.processType = root.PandaAnalyzer.kTT
    fin = root.TFile.Open(fullPath)

    print fullPath
    print fin

    tree = fin.FindObjectAny("events")
    infotree = fin.FindObjectAny("all")
    print tree,infotree

    skimmer.SetDataDir(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/')
    skimmer.SetOutputFile('testskim.root')
    skimmer.Init(tree,infotree)

    skimmer.Run()
    print 'done running'
    skimmer.Terminate()
    print 'done terminating'

  fn(argv[1]) 
