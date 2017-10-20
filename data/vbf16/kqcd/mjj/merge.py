#!/usr/bin/env python

# merge efficiencies from https://gitlab.cern.ch/cms-vbf-invisible/mcweights/tree/master/kfactor_vjet_qcd/latest 

from sys import argv
infile = argv[1]
outfile = argv[2]
base = argv[1]
output = argv[2]

import ROOT as root
from array import array 

f_out = root.TFile(output, 'UPDATE')

directory = 'kfactors_shape'
bins = [200, 500, 1000, 1500, 5000]
bins_mjj = array('f', bins)

f_in = root.TFile(base)
f_in.cd(directory)

h_merged = None

for ib in xrange(len(bins) - 1):
    r = (bins[ib], bins[ib+1])
    h1 = f_in.Get(directory + '/kfactor_vbf_mjj_%i_%i_smoothed' % r)
    if not h_merged:
        bins_pt = array('f', [h1.GetBinLowEdge(i) for i in xrange(1, h1.GetNbinsX()+2)])
        h_merged = root.TH2F('merged', 'merged', len(bins_pt)-1, bins_pt, len(bins_mjj)-1, bins_mjj)
    for jb in xrange(1, h1.GetNbinsX()+1):
        h_merged.SetBinContent(jb, ib + 1, h1.GetBinContent(jb))

f_out.WriteTObject(h_merged, 'h_%s'%directory, 'overwrite')

directory = 'kfactors_cc'
h1 = f_in.Get(directory + '/kfactor_vbf')
f_out.WriteTObject(h1, 'h_%s'%directory, 'overwrite')

f_out.Close()
