#!/usr/bin/env python

from sys import argv
from os import getenv
which = argv[1]
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
Load('Tools','BranchAdder')

ba = root.BranchAdder()
ba.formula = 'max(0,min(40,npv))'
ba.newBranchName = 'sf_pu2016_fixed'

fin = root.TFile(which,'UPDATE')
tin = fin.Get('events')
fpu = root.TFile('/data/t3home000/bmaier/flat_v9/cr_dimuon/puWeight.root')
hpu = fpu.Get('hPU')


ba.AddBranchFromHistogram(tin,hpu)
fin.WriteTObject(tin,'events','Overwrite')
fin.Close()
