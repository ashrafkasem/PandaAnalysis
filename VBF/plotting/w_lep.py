#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange

basedir = getenv('PANDA_FLATDIR') + '/../vbf_004_6w/'
infile = basedir+'/WJets.root'


parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
args = parser.parse_args()

figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Tools.Load import Load
from PandaCore.Tools.Misc import *
from PandaCore.Tools.root_interface import Selector
Load('HistogramDrawer')

lumi=35900
import PandaAnalysis.VBF.PandaSelection as sel 
base_cut = tAND(sel.cnc, sel.cuts['signal'])
weight = sel.weights['signal']%lumi

plot = root.HistogramDrawer()
plot.SetRatioStyle()
plot.SetLumi(lumi/1000.)
plot.DrawEmpty(True)
plot.AddCMSLabel()
plot.AddLumiLabel()
#plot.SetAutoRange(False)
plot.InitLegend()
root.gStyle.SetOptStat(0)


s = Selector()
f = root.TFile.Open(infile); t = f.Get('events')
branches = ['pfmet', 'jot12Mass', '1',
            'fabs(genTauEta)', 'fabs(genMuonEta)', 'fabs(genElectronEta)', 
            'genTauPt', 'genElectronPt', 'genMuonPt',
            weight]
s.read_tree(t, branches = branches, cut = base_cut)

t_mask = ((s['genTauPt'] > 1) &
          (s['genMuonPt'] < 5) &
          (s['genElectronPt'] < 5))
e_mask = ( #(s['genMuonPt'] < 5) &
          (s['genElectronPt'] > 5))
m_mask = ((s['genMuonPt'] > 5) &
          (s['genElectronPt'] < 5))


ratios_to_plot = {
  'pfmet' : ([250, 350, 500, 700, 1000],'E_{T}^{miss} [GeV]'),
  '1' : ([0, 2],'yield'),

  'jot12Mass' : ([200, 600, 1000, 1500, 2500, 3500, 5000],'m_{jj} [GeV]')
}
for k,v in ratios_to_plot.iteritems():
    h_inc = s.draw(k, weight, vbins=v[0])

    h_t = s.draw(k, weight, mask=t_mask, vbins=v[0])
    h_e = s.draw(k, weight, mask=e_mask, vbins=v[0])
    h_m = s.draw(k, weight, mask=m_mask, vbins=v[0])

    plot.Reset(False)
    for h in [h_t, h_e, h_m]:
        h.Divide(h_inc)
        h.GetXaxis().SetTitle(v[1])
        h.GetYaxis().SetTitle('Fraction of W events')
        h.SetMaximum(0.8)

    plot.AddHistogram(h_t, 'W#rightarrow#tau', 5, 1, 'hist')
    plot.AddHistogram(h_m, 'W#rightarrow#mu', 6, 2, 'hist')
    plot.AddHistogram(h_e, 'W#rightarrowe', 7, 3, 'hist')

    plot.Draw(args.outdir,'ratio_mjj_'+k)


# now do lepton etas
plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.SetLumi(lumi/1000.)
plot.DrawEmpty(True)
plot.AddCMSLabel()
plot.AddLumiLabel()
plot.SetAutoRange(True)
plot.InitLegend()
root.gStyle.SetOptStat(0)

eta_bins = (10, 0, 5)
h_t = s.draw('fabs(genTauEta)', weight, mask=t_mask, fbins=eta_bins)
h_m = s.draw('fabs(genMuonEta)', weight, mask=m_mask, fbins=eta_bins)
h_e = s.draw('fabs(genElectronEta)', weight, mask=e_mask, fbins=eta_bins)

plot.AddHistogram(h_t, 'W#rightarrow#tau', 5, 1, 'hist')
plot.AddHistogram(h_m, 'W#rightarrow#mu', 6, 2, 'hist')
plot.AddHistogram(h_e, 'W#rightarrowe', 7, 3, 'hist')
plot.Draw(args.outdir,'eta_mjj')

print '%3s fraction beyond accp = %.2f%%'%('tau', 100*h_t.Integral(h_t.FindBin(2.5),20)/h_t.Integral())
print '%3s fraction beyond accp = %.2f%%'%('mu', 100*h_m.Integral(h_m.FindBin(2.5),20)/h_m.Integral())
print '%3s fraction beyond accp = %.2f%%'%('e', 100*h_e.Integral(h_e.FindBin(2.5),20)/h_e.Integral())
