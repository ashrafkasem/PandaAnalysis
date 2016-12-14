#!/usr/bin/env python

from re import sub
from sys import argv,exit
inpath = argv[1]
argv=[]

from os import system,getenv
from PandaCore.Tools.Load import *
import ROOT as root

if __name__ == "__main__":

  Load('PandaAnalysisFlat','GenAnalyzer')
  
  def fn(fullPath):

    skimmer = root.GenAnalyzer()
   
    skimmer.firstEvent=0
    skimmer.lastEvent=10000
    fin = root.TFile.Open(fullPath)

    print fullPath

    tree = fin.FindObjectAny("genEvents")

    skimmer.SetOutputPath('./')
    skimmer.order = root.GenAnalyzer.kNLO
    skimmer.processType = root.GenAnalyzer.kW
    skimmer.AddInput(tree,1.5,'w_nlo')

    skimmer.Run()
    print 'done running'
    skimmer.Terminate()
    print 'done terminating'

  fn(inpath) 
