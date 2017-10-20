#!/usr/bin/env python

from sys import argv
from os import getenv
which = argv[1]
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
Load('BranchAdder')

ba = root.BranchAdder()
ba.formula = '1'
ba.newBranchName = 'sf_qcdV_VBFTight'

fin = root.TFile(getenv('PANDA_FLATDIR')+'/'+which,'UPDATE')
tin = fin.Get('events')
ba.AddBranchFromFormula(tin)
fin.WriteTObject(tin,'events','Overwrite')
fin.Close()
