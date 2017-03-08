#!/usr/bin/env python

from sys import argv 
import argparse

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--label',metavar='label',type=str)
args = parser.parse_args()
argv = []

import ROOT as root
from PandaCore.Tools.Load import Load

Load('Tools','EventSyncher')
Load('Drawers','CanvasDrawer')

es = root.EventSyncher()

if args.label=='MET':
	f1 = '/mnt/hadoop/scratch/bmaier/panda/v7/flat/MET.root'
	f2 = '/mnt/hadoop/scratch/bmaier/panda/v9/flat/MET.root'
else:
	f1 = '/mnt/hadoop/scratch/bmaier/panda/v7/flat/TTbar_CUETP8M2T4.root'
	f2 = '/mnt/hadoop/scratch/bmaier/panda/v9/flat/TTbar.root'
es.preselection = '(( metFilter==1 ) && ( (( nFatjet==1 && fj1Pt>250 && fabs(fj1Eta)<2.4 ) && ( pfUWmag>250 && dphipfUW>0.5 && nLoosePhoton==0 && nTau==0 && nLooseLep==1 && looseLep1IsTight==1 && abs(looseLep1PdgId)==13 )) ))' 

formulae = [
		('fj1MSD',0,300,'m_{SD}'),
		('fj1Pt',250,500,'p_{T}'),
		('fj1RawPt',250,500,'raw p_{T}'),
		]

for f in formulae:
	es.AddFormula(*f)

es.RunFiles(f1,f2,'events')

hists = es.PlayViolin(0.2)

plot = root.CanvasDrawer()
plot.SetTDRStyle()
c = root.TCanvas()

for iF in xrange(len(formulae)):
	hist = hists[iF]
	c.Clear()
	hist.Draw('VIOLIN')
	plot.AddCMSLabel(.16,.94)
	plot.SetLumi(36.6); plot.AddLumiLabel(True)
	plot.SetCanvas(c)
	plot.Draw(args.outdir+'/',args.label + '_' + formulae[iF][0])
