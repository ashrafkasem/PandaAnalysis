#!/usr/bin/env python

# merge efficiencies from https://gitlab.cern.ch/cms-vbf-invisible/mcweights/tree/master/met_trigger_weight

from sys import argv
base = argv[1]
output = argv[2]

import ROOT as root
from numpy import arange 
from array import array 

bins = [0, 800, 1200, 1700, 3000]

f_out = root.TFile(output, 'RECREATE')
h_merged = None

for ib in xrange(len(bins) - 1):
    r = (bins[ib], bins[ib+1])
    f = root.TFile.Open(base + '%i.0_%i.0.root' % r)
    eff = f.Get('efficiency')
    eff_fn = f.Get('efficiency_func')
    nbx = eff.GetPassedHistogram().GetNbinsX()
    if not h_merged:
        hdummy= eff.GetPassedHistogram()
        x_arr = array('f', arange(hdummy.GetBinLowEdge(1), hdummy.GetBinLowEdge(nbx+1), 1))
        # x_arr = array('f', [.GetXaxis().GetXbins().At(i) for i in xrange(nbx)])
        y_arr = array('f', bins)
        x_n = len(x_arr) - 1
        y_n = len(y_arr) - 1
        h_merged = root.TH2F('h', 'h', 
                             x_n, x_arr,
                             y_n, y_arr)
        h_merged.SetDirectory(0)
    xaxis = h_merged.GetXaxis()        
    for ibx in xrange(1, x_n + 1):
        x_val = xaxis.GetBinCenter(ibx)
        # h_merged.SetBinContent(ibx, ib + 1, eff.GetEfficiency(ibx))
        h_merged.SetBinContent(ibx, ib + 1, eff_fn.Eval(x_val))
    f_out.WriteTObject(eff, 'eff_%i_%i' % r)
    f_out.WriteTObject(eff_fn, 'fn_%i_%i' % r)

f_out.WriteTObject(h_merged, 'h_eff')
f_out.Close()
