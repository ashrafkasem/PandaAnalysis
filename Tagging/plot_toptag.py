#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
import argparse

basedir = getenv('PANDA_FLATDIR')
figsdir = basedir+'/figs'

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--indir',metavar='indir',type=str,default=basedir)
parser.add_argument('--outdir',metavar='outdir',type=str,default=figsdir)
parser.add_argument('--cut',metavar='cut',type=str,default=None)
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
cut = 'nFatjet==1 && fj1Pt>250 && fj1MaxCSV>0.46 && nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && UWmag>250 && isojetNBtags>0'
if not args.cut:
  label = 'noCut_'
  plotlabel = None
elif args.cut=='mass':
  cut = tAND(cut,'fj1MSD>110 && fj1MSD<210')
  label = 'massCut_'
  plotlabel = '110 < m_{SD} < 210 GeV'

### LOAD PLOTTING UTILITY ###
plot = root.PlotUtility()
plot.Stack(True)
plot.InitLegend()
plot.SetCut(tcut(cut))
plot.Ratio(1) 
plot.SetLumi(lumi/1000)
plot.DrawMCErrors(True)
plot.SetTDRStyle()
plot.SetNormFactor(False)
plot.AddCMSLabel()
plot.AddLumiLabel()
if plotlabel:
  plot.AddPlotLabel(plotlabel,.18,.77,False,42,.04)

weight = '0.1*%f*normalizedWeight*sf_pu*sf_lep*sf_ewkV*sf_qcdV*sf_sjbtag1*sf_btag1*sf_tt'%lumi
plot.SetMCWeight(weight)

### DEFINE PROCESSES ###
wjets     = root.Process('W+jets',root.kWjets)
diboson   = root.Process('Diboson',root.kDiboson)
ttbar     = root.Process('t#bar{t} [matched]',root.kTTbar); ttbar.additionalCut = root.TCut('(fj1IsMatched==1&&fj1GenSize<1.2)')
ttbarunmatched     = root.Process('t#bar{t} [unmatched]',root.kExtra1); ttbarunmatched.additionalCut = root.TCut('(fj1IsMatched==0||fj1GenSize>1.2)')
singletop = root.Process('Single t',root.kST)
qcd       = root.Process("QCD",root.kQCD)
data      = root.Process("Data",root.kData); data.additionalCut = root.TCut('(trigger&1)!=0')
processes = [diboson,singletop,wjets,ttbarunmatched,ttbar]

### ASSIGN FILES TO PROCESSES ###
wjets.AddFile(basedir+'WJets.root')
diboson.AddFile(basedir+'Diboson.root')
ttbar.AddFile(basedir+'TTbar.root')
ttbarunmatched.AddFile(basedir+'TTbar.root')
singletop.AddFile(basedir+'SingleTop.root')
qcd.AddFile(basedir+'QCD.root')

data.AddFile(basedir+'MET.root') 
processes.append(data)

for p in processes:
  plot.AddProcess(p)

plot.AddDistribution(root.Distribution('fj1ECFN_2_4_20/pow(fj1ECFN_1_3_20,2)',0,5,20,'N_{3}(#beta=2.0)','Events',999,-999,'N3_20'))

plot.AddDistribution(root.Distribution('fj1ECFN_2_4_10/pow(fj1ECFN_1_3_10,2)',0.5,3.5,20,'N_{3}(#beta=1.0)','Events',999,-999,'N3_10'))

plot.AddDistribution(root.Distribution('fj1ECFN_2_4_40/pow(fj1ECFN_1_3_40,2)',0,5,20,'N_{3}(#beta=4.0)','Events',999,-999,'N3_40'))

plot.AddDistribution(root.Distribution('fj1MSD',50,300,20,'fatjet m_{SD} [GeV]','Events/12.5 GeV'))

### DRAW AND CATALOGUE ###
plot.DrawAll(figsdir+'/'+label)
