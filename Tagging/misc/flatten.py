#!/usr/bin/env python

from sys import argv
from os import getenv
which = argv[1]
argv = []

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *

Load('BranchAdder')

def gethisto(tree,formula,binlo,binhi,additionalcut=None):
  nbins=50
  h = root.TH1D('h','h',nbins,binlo,binhi)
  if additionalcut:
    tree.Draw(formula+'>>h',tTIMES('normalizedWeight',additionalcut))
  else:
    tree.Draw(formula+'>>h','normalizedWeight')

  for iB in xrange(1,nbins+1):
    val = h.GetBinContent(iB);
    if val:
      h.SetBinContent(iB,1./val)
    else:
      print formula,h.GetBinCenter(iB)

  return h


def addbranchesFormula(fpath,additionalcut=None):
  fin = root.TFile(fpath,'UPDATE')
  jets = fin.Get('events')
  ba = root.BranchAdder()

  ba.formula = '1./(250*TMath::Exp((232-fj1Pt)*0.0396)+200*TMath::Exp((235-fj1Pt)*0.0157)+TMath::Exp((583.-fj1Pt)*0.00672))'
  ba.newBranchName = 'ptweight_qcd'
  ba.AddBranchFromFormula(jets)

  fin.WriteTObject(jets,'events','Overwrite')
  fin.Close()

def addbranches(fpath,additionalcut=None):
  fin = root.TFile(fpath,'UPDATE')
  jets = fin.Get('events')
  ba = root.BranchAdder()

  hpt = gethisto(jets,'fj1Pt',250,1000,additionalcut)
  ba.formula = 'fj1Pt'
  ba.newBranchName = 'ptweight_binned'
  ba.AddBranchFromHistogram(jets,hpt)

  fin.WriteTObject(jets,'events','Overwrite')
  fin.Close()

basedir = getenv('PANDA_FLATDIR')

if which=='ZpTT':
  addbranches(basedir+'/ZpTT.root','fj1IsMatched==1&&fj1GenSize<1.44')
elif which=='ZpWW':
  addbranches(basedir+'/ZpWW.root','fj1IsMatched==1&&fj1GenSize<1.44')
elif which=="ZpA0h":
  addbranches(basedir+'/ZpA0h.root','fj1IsMatched==1&&fj1GenSize<1.44')
else:
  addbranchesFormula(basedir+'/'+which+'.root')
