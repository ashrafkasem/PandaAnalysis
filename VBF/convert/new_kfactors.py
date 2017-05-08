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
h_num = fk.Get('bosonPt_NLO_vbf')
h_den = fk.Get('bosonPt_LO_vbf')
h_ratio = h_num.Clone(); h_ratio.Divide(h_den)

ba.formula = 'genBosonPt'
ba.newBranchName = 'sf_qcdVVBF'

fin = root.TFile(getenv('PANDA_FLATDIR')+which+'.root','UPDATE')
tin = fin.Get('events')
ba.AddBranchFromHistogram(tin,h_ratio)
fin.WriteTObject(tin,'events','Overwrite')
fin.Close()
