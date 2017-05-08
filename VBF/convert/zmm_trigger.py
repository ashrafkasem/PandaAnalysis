#!/usr/bin/env python

from sys import argv
from os import getenv
which = argv[1]
trigger = argv[2]
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
Load('BranchAdder')

ba = root.BranchAdder()

fk = root.TFile(trigger)
h_ratio = fk.Get('hden_monojet_recoil_clone_passed') 

ba.formula = 'pfUZmag'
ba.newBranchName = 'sf_metTrigZmm'

fin = root.TFile(getenv('PANDA_FLATDIR')+which+'.root','UPDATE')
tin = fin.Get('events')
ba.AddBranchFromHistogram(tin,h_ratio)
fin.WriteTObject(tin,'events','Overwrite')
fin.Close()
