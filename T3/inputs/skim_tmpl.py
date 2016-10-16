#!/usr/bin/env python

from re import sub
from sys import argv,exit
from os import system,getenv,path
from time import clock
import subprocess

which = int(argv[1])
sname = argv[0]
argv=[]

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import *

if __name__ == "__main__":
  
  Load('PandaAnalysisFlat','PandaAnalyzer')

  def fn(shortName,longName,isData,xsec):
    outfilename = shortName+'.root'

    start=clock()

    eosPath = 'root://eoscms.cern.ch//store/user/snarayan'
    eosEXOPath = 'root://eoscms.cern.ch//store/group/phys_exotica'
    cernboxPath = 'root://eosuser//eos/user/s/snarayan'
    cernboxBPath = 'root://eosuser//eos/user/b/bmaier'
    fullPath = sub(r'\${CERNBOXB}',cernboxBPath,
        sub(r'\${CERNBOX}',cernboxPath,
          sub(r'\${EOS}',eosPath,
            sub(r'\${EOS2}',eosEXOPath,longName))))
    PInfo(sname,fullPath)

    system('xrdcp %s input.root'%fullPath)

    skimmer = root.PandaAnalyzer()
    skimmer.isData=isData
    skimmer.applyJson=False
    skimmer.SetPreselectionBit(root.PandaAnalyzer.kMonotop)
    processType=root.PandaAnalyzer.kNone
    if not isData:
      if 'ZJets' in fullPath or 'DY' in fullPath:
        processType=root.PandaAnalyzer.kZ
      elif 'WJets' in fullPath:
        processType=root.PandaAnalyzer.kW
      elif 'GJets' in fullPath:
        processType=root.PandaAnalyzer.kA
      elif 'TTJets' in fullPath or 'TT_' in fullPath:
        processType=root.PandaAnalyzer.kTT
    skimmer.processType=processType 

    try:
      fin = root.TFile.Open('input.root')
      tree = fin.FindObjectAny("events")
      infotree = fin.FindObjectAny("all")
    except:
      exit(2) # file open error => xrootd?

    skimmer.SetDataDir(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/')
    skimmer.SetOutputFile('output.root')
    skimmer.Init(tree,infotree)

    skimmer.Run()
    skimmer.Terminate()

    mvargs = '$PWD/output.root XXXX%s'%outfilename
    system(mvargs)
    system('rm input.root')

    PInfo(sname,'finished in %f'%(clock()-start)); start=clock()

  
  cfg = open('local.cfg')
  lines = list(cfg)
  ll = lines[which].split()
  shortname = ll[0]
  isData = (ll[1]!="MC")
  xsec = float(ll[2])
  longname = ll[3]
  cfg.close() # free up memory
  del lines
  fn(shortname,longname,isData,xsec)

