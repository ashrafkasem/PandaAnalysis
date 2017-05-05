#!/usr/bin/env python

import sys
import glob

pattern = sys.argv[1]

sys.argv = []

import ROOT as root

files = [pattern]

for fpath in files:
    f = root.TFile.Open(fpath)
    t = f.Get('events')
    n_rw = sum([('rw_' in b.GetName()) for b in t.GetListOfBranches()])
#    print n_rw
    if n_rw!=377 and n_rw!=22:
        print n_rw, fpath
