#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv,exit
import argparse

basedir = getenv('PANDA_FLATDIR')
figsdir = basedir+'/figs'

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--indir',metavar='indir',type=str,default=basedir)
parser.add_argument('--outdir',metavar='outdir',type=str,default=figsdir)
parser.add_argument('--sel',metavar='sel',type=str,default='pass')
args = parser.parse_args()

figsdir = args.outdir
basedir = args.indir
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
tcut = root.TCut
Load('Drawers','PlotUtility')

### SET GLOBAL VARIABLES ###
lumi = 12918.
logy=False
nlo = 'sf_ewkV*sf_qcdV'
cut = 'nFatjet==1 && fj1Pt>250 && fj1MaxCSV>0.46 && nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && UWmag>250 && isojetNBtags==1'
weight = '%f*normalizedWeight*sf_pu*sf_lep*%s*sf_sjbtag1*sf_btag1*sf_tt*sf_metTrig'%(lumi,nlo)
if args.sel=='pass':
  label = 'pass_'
  cut = tAND(cut,'top_ecfv8_bdt>0.66')
else:
  label = 'fail_'
  cut = tAND(cut,'top_ecfv8_bdt<0.66')

### LOAD PLOTTING UTILITY ###
plot = root.PlotUtility()
plot.Stack(True)
plot.InitLegend()
plot.SetCut(tcut(cut))
plot.Ratio(1) 
plot.FixRatio(.5)
plot.SetLumi(lumi/1000)
plot.DrawMCErrors(True)
plot.SetTDRStyle()
plot.SetNormFactor(False)
plot.AddCMSLabel()
plot.AddLumiLabel()
plot.SetMCWeight(weight)


### DEFINE PROCESSES ###
wjets     = root.Process('W',root.kWjets)
others   = root.Process('others',root.kDiboson)
ttbar    = root.Process('tt_matched',root.kTTbar); ttbar.additionalCut = root.TCut('(fj1IsMatched==1&&fj1GenSize<1.44)')
ttbarunmatched     = root.Process('tt_unmatched',root.kExtra1); ttbarunmatched.additionalCut = root.TCut('(fj1IsMatched==0||fj1GenSize>1.44)')
data      = root.Process("Data",root.kData)
data.additionalCut = root.TCut('(trigger&1)!=0')
processes = [others,wjets,ttbarunmatched,ttbar]

### ASSIGN FILES TO PROCESSES ###
data.AddFile(basedir+'MET.root') 
# others.AddFile(basedir+'QCD.root')
wjets.AddFile(basedir+'WJets.root')
others.AddFile(basedir+'Diboson.root')
ttbar.AddFile(basedir+'TTbar.root')
ttbarunmatched.AddFile(basedir+'TTbar.root')
others.AddFile(basedir+'SingleTop.root')
processes.append(data)

for p in processes:
  plot.AddProcess(p)


plot.AddDistribution(root.Distribution('fj1MSD',40,450,30,'fatjet m_{SD} [GeV]','Events/14 GeV'))

plot.DrawAll(figsdir+'/'+label)
