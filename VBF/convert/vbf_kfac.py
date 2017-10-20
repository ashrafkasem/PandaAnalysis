#!/usr/bin/env python

from sys import argv
from os import getenv
which = argv[1]
kfactor = argv[2]
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
Load('BranchAdder')

ba = root.BranchAdder()

fk = root.TFile(kfactor)
h = fk.Get('TH2F_kFactor')

if 'WJets' in which:
    ba.formula = 'pfUAmag'
elif 'ZJets' in which:
    ba.formula = 'pfUZmag'
else:
    ba.formula = 'pfmet'
ba.formulaY = 'jot12Mass'
ba.newBranchName = 'sf_qcdVVBF_v2'

fin = root.TFile(getenv('PANDA_FLATDIR')+which+'.root','UPDATE')
tin = fin.Get('events')
ba.AddBranchFromHistogram2D(tin,h)
fin.WriteTObject(tin,'events','Overwrite')
fin.Close()
