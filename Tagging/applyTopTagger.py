#!/usr/bin/env python

from os import getenv
from sys import argv

import cfgModeled as cfg
#import cfgTopModeled as cfg

sample = argv[1]
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *

Load('Learning','TMVABranchAdder')

workdir = getenv('PANDA_FLATDIR')

ba = root.TMVABranchAdder()
ba.treename='events'
ba.defaultValue=-1.2
ba.presel='fj1ECFN_2_4_20>0'

for v in cfg.variables:
  ba.AddVariable(v[0],v[2])

for v in cfg.formulae:
  ba.AddFormula(v[0],v[2])

for s in cfg.spectators:
  ba.AddSpectator(s[0])

ba.BookMVA('top_ecfv6fixed_bdt',workdir+'/training/top_ecfbdt_v6_BDT.weights.xml')

ba.RunFile(workdir+'/'+sample+'.root')

