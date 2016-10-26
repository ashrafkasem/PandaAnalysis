#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv,exit
import argparse

basedir = getenv('PANDA_FLATDIR')
figsdir = basedir+'/figs'

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--indir',metavar='indir',type=str,default=basedir)
parser.add_argument('--outdir',metavar='outdir',type=str,default=figsdir)
parser.add_argument('--disc',metavar='disc',type=str,default='top_ecfv7_bdt')
parser.add_argument('--tagged',metavar='tagged',type=str,default='True')
parser.add_argument('--sel',metavar='sel',type=str,default='tag')
parser.add_argument('--cut',metavar='cut',type=float,default=0.5)
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
if args.sel=='mistag':
  cut = 'nFatjet==1 && fj1Pt>250 && fj1MaxCSV<0.46 && nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && UWmag>250 && isojetNBtags==0'
  weight = '%f*normalizedWeight*sf_pu*sf_lep*%s*sf_sjbtag0*sf_btag0*sf_tt*sf_metTrig'%(lumi,nlo)
  label = 'mistag_'
elif args.sel=='photon':
  cut = 'nFatjet==1 && fj1Pt>250 && nLooseLep==0 && nLoosePhoton==1 && loosePho1IsTight==1 && nTau==0 && UAmag>250'
  weight = '%f*normalizedWeight*sf_pu*sf_lep*%s*sf_tt*sf_phoTrig'%(lumi,nlo)
  label = 'photon_'
else:
  cut = 'nFatjet==1 && fj1Pt>250 && fj1MaxCSV>0.46 && nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && UWmag>250 && isojetNBtags==1'
  weight = '%f*normalizedWeight*sf_pu*sf_lep*%s*sf_sjbtag1*sf_btag1*sf_tt*sf_metTrig'%(lumi,nlo)
  label = 'tag_'

disclabels = {
    'fj1Tau32SD':'#tau_{32}^{SD}',
    'top_ecfv6_bdt':'ECF BDT',
    'top_ecfv7_bdt':'ECF+#tau_{32}^{SD} BDT',
    'top_ecfv8_bdt':'ECF+#tau_{32}^{SD}+f_{rec} BDT',
    }

label += args.disc+'_'
if args.tagged.upper()=='TRUE':
  label += 'pass_'
  if 'Tau' in args.disc:
    cut = tAND(cut,'%s<%f'%(args.disc,args.cut))
    plotlabel = '%s<%.2f'%(disclabels[args.disc],args.cut)
  else:
    cut = tAND(cut,'%s>%f'%(args.disc,args.cut))
    plotlabel = '%s>%.2f'%(disclabels[args.disc],args.cut)
else:
  label += 'fail_'
  if 'Tau' in args.disc:
    cut = tAND(cut,'%s>%f'%(args.disc,args.cut))
    plotlabel = '%s>%.2f'%(disclabels[args.disc],args.cut)
  else:
    cut = tAND(cut,'%s<%f'%(args.disc,args.cut))
    plotlabel = '%s<%.2f'%(disclabels[args.disc],args.cut)

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
if plotlabel:
  plot.AddPlotLabel(plotlabel,.18,.77,False,42,.04)
plot.SetMCWeight(weight)

### DEFINE PROCESSES ###
wjetsq     = root.Process('W+q',root.kWjets); wjetsq.additionalCut = root.TCut('abs(fj1HighestPtGen)!=21')
wjetsg     = root.Process('W+g',root.kExtra2); wjetsg.additionalCut = root.TCut('abs(fj1HighestPtGen)==21')
diboson   = root.Process('Diboson',root.kDiboson)
ttbar     = root.Process('t#bar{t} [matched]',root.kTTbar); ttbar.additionalCut = root.TCut('(fj1IsMatched==1&&fj1GenSize<1.44)')
ttbarunmatched     = root.Process('t#bar{t} [unmatched]',root.kExtra1); ttbarunmatched.additionalCut = root.TCut('(fj1IsMatched==0||fj1GenSize>1.44)')
singletop = root.Process('Single t',root.kST)
qcd       = root.Process("QCD",root.kQCD)
data      = root.Process("Data",root.kData)
gjetsq    = root.Process('#gamma+q',root.kGjets); gjetsq.additionalCut = root.TCut('abs(fj1HighestPtGen)!=21')
gjetsg    = root.Process('#gamma+g',root.kExtra3); gjetsg.additionalCut = root.TCut('abs(fj1HighestPtGen)==21')
if args.sel=='photon':
  data.additionalCut = root.TCut('(trigger&4)!=0')
  qcd.useCommonWeight=False
  qcd.additionalWeight = root.TCut('sf_phoPurity')
  qcd.additionalCut = root.TCut('(trigger&4)!=0')
  processes = [qcd,gjetsg,gjetsq]
else:
  data.additionalCut = root.TCut('(trigger&1)!=0')
  processes = [diboson,singletop,wjetsg,wjetsq,ttbarunmatched,ttbar]

### ASSIGN FILES TO PROCESSES ###
if args.sel=='photon':
  gjetsq.AddFile(basedir+'GJets.root')
  gjetsg.AddFile(basedir+'GJets.root')
  data.AddFile(basedir+'SinglePhoton.root') 
  qcd.AddFile(basedir+'SinglePhoton.root') 
else:
  data.AddFile(basedir+'MET.root') 
  qcd.AddFile(basedir+'QCD.root')
  wjetsq.AddFile(basedir+'WJets.root')
  wjetsg.AddFile(basedir+'WJets.root')
  diboson.AddFile(basedir+'Diboson.root')
  ttbar.AddFile(basedir+'TTbar.root')
  ttbarunmatched.AddFile(basedir+'TTbar.root')
  singletop.AddFile(basedir+'SingleTop.root')
processes.append(data)

for p in processes:
  plot.AddProcess(p)

plot.AddDistribution(root.Distribution('fj1MSD',40,450,20,'fatjet m_{SD} [GeV]','Events/12.5 GeV'))

### DRAW AND CATALOGUE ###
plot.DrawAll(figsdir+'/'+label)
