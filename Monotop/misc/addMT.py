#!/usr/bin/env python

from sys import argv
from os import getenv
which = argv[1]
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
Load('Tools','BranchAdder')

ba = root.BranchAdder()
ba.formula = 'MT(looseLep1Pt,looseLep1Phi,pfmet,pfmetphi)'
ba.newBranchName = 'fixed_mt'

fin = root.TFile(getenv('PANDA_FLATDIR')+'/'+which+'.root','UPDATE')
tin = fin.Get('events')
ba.AddBranchFromFormula(tin)
fin.WriteTObject(tin,'events','Overwrite')
fin.Close()
