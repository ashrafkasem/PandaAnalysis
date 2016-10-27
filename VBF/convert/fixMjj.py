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
ba.formula = 'Mxx(jot1Pt,jot1Eta,jot1Phi,jot1M,jot2Pt,jot2Eta,jot2Phi,jot2M)'
ba.newBranchName = 'fixed_mjj'

fin = root.TFile(getenv('PANDA_ZEYNEPDIR')+'/merged/'+which+'.root')
tin = fin.Get('events')
ba.AddBranchFromFormula(tin)
fin.WriteTObject(tin,'events','Overwrite')
fin.Close()
