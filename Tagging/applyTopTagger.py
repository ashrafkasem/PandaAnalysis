#!/usr/bin/env python

from os import getenv
from sys import argv

#import cfgModeled as cfg
import cfgAll as cfg

nextra = int(argv[1])
sample = argv[2]
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

for v in cfg.variables[:nextra]:
  ba.AddVariable(v[0],v[2])

for v in cfg.formulae[:12]:
  ba.AddFormula(v[0],v[2])

for s in cfg.spectators:
  ba.AddSpectator(s[0])

#ba.BookMVA('top_ecfv6_bdt',workdir+'/training/top_ecfbdt_v6_BDT.weights.xml')
#ba.BookMVA('top_ecfv7_bdt',workdir+'/training/top_ecfbdt_v7_BDT.weights.xml')
#ba.BookMVA('top_ecfv8_bdt',workdir+'/training/top_ecfbdt_v8_BDT.weights.xml')
v = nextra+12
ba.BookMVA('top_ecfv%i_bdt'%v,workdir+'/training/top_ecfbdt_v%i_BDT.weights.xml'%v)

ba.RunFile(workdir+'/'+sample+'.root')

