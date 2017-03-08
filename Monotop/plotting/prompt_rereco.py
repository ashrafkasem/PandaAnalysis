#!/usr/bin/env python

from os import system,getenv
from sys import argv,exit
import argparse
from collections import namedtuple

### SET GLOBAL VARIABLES ###
rerecoDir = '/home/snarayan/home000/store/panda/v_8020_2_3/' 
promptDir = '/home/snarayan/home000/store/panda/v_8020_2_PR/'
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
parser.add_argument('--region',metavar='region',type=str)
parser.add_argument('--era',metavar='era',type=str)
args = parser.parse_args()
lumi = 35000

region = args.region
era = args.era
sname = argv[0]
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.Monotop.NoMassPFSelection as sel
Load('Drawers','PlotUtility')

Era = namedtuple('Era',['lo','hi','lumi'])

all_eras = {
    'B' : Era(272007,275376,5927.9),
    'C' : Era(275657,276283,2628.9),
    'D' : Era(276315,276811,4341.9),
    'E' : Era(276831,277420,4040.3),
    'F' : Era(277772,278808,3160.1),
    'G' : Era(278820,280385,7552.5),
    'H' : Era(280919,284044,8761.8),
    }

run_boundaries = []
for e in era:
  lo,hi = all_runs[e]
  run_boundaries += [lo,hi]
runs = (min(run_boundaries),max(run_boundaries))
cut = tAND(tAND(sel.cuts[args.region],args.cut),'runNumber>%i && runNumber<%i'%(runs[0],runs[1]))

### LOAD PLOTTING UTILITY ###
plot = root.PlotUtility()
plot.Stack(True)
plot.DrawEmpty(True)
plot.SetNormFactor(True)
plot.Ratio(True)
plot.FixRatio(1)
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
#plot.AddLumiLabel(True)
plot.SetEvtNum("eventNumber")
if region=='signal':
  cut = tAND(cut,'eventNumber%5==0')
  plot.AddPlotLabel('Every 5th event',.18,.7,False,42,.04)
plot.SetCut(root.TCut(cut))
plot.SetMCWeight(root.TCut('1'))

region_labels = {
    'signal' : 'Signal region',
    'wen' : 'Single electron CR',
    'zee' : 'Dielectron CR',
    'wmn' : 'Single muon CR',
    'zmm' : 'Dimuon CR',
    }

plot.AddPlotLabel('#splitline{%s}{Run2016%s}'%(region_labels[region],era),.18,.77,False,42,.04)

### DEFINE PROCESSES ###
rereco      = root.Process("ReReco",root.kData)
prompt      = root.Process("PromptReco",root.kExtra1)
processes = [prompt,rereco]

### ASSIGN FILES TO PROCESSES ###
if 'signal' in region or 'muon' in region:
  rereco.AddFile(rerecoDir+'MET.root')
  prompt.AddFile(promptDir+'MET.root')
  lep='#mu'
elif 'electron' in region:
  rereco.AddFile(rerecoDir+'SingleElectron.root')
  prompt.AddFile(promptDir+'SingleElectron.root')
  lep='e'

for p in processes:
  plot.AddProcess(p)

#recoilBins = [200,250,300,350,400,500,600,750,1000]
#recoilBins = [200,250,300,350,400,500,600,1000]
recoilBins = [200., 230., 260.0, 290.0, 320.0, 350.0, 390.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]
nRecoilBins = len(recoilBins)-1

### CHOOSE DISTRIBUTIONS, LABELS ###
plot.AddDistribution(root.Distribution('pfmetphi',-3.142,3.142,20,'PF MET #phi','Events'))
#plot.AddDistribution(root.Distribution('top_ecf_bdt',-1,1,20,'ECF+#tau_{32}^{SD}+f_{rec} BDT','Events'))
plot.AddDistribution(root.Distribution('fj1Tau32SD',0,1,20,'Groomed #tau_{32}','Events',999,-999,'tau32SD'))
plot.AddDistribution(root.Distribution('fj1Tau32',0,1,20,'#tau_{32}','Events',999,-999,'tau32'))
plot.AddDistribution(root.Distribution('jet1Pt',15,1000,20,'leading jet p_{T} [GeV]','Events'))
# plot.AddDistribution(root.Distribution('nJet',-0.5,8.5,9,'N_{jet}','Events'))
plot.AddDistribution(root.Distribution('fj1Pt',250,1000,20,'fatjet p_{T} [GeV]','Events/37.5 GeV'))
plot.AddDistribution(root.Distribution('fj1Eta',-2.5,2.5,20,'fatjet #eta','Events'))
plot.AddDistribution(root.Distribution('fj1MSD',50,250,20,'fatjet m_{SD} [GeV]','Events/10 GeV'))
#plot.AddDistribution(root.Distribution('fj1MSD',0,450,20,'fatjet m_{SD} [GeV]','Events/12.5 GeV'))
plot.AddDistribution(root.Distribution('fj1MaxCSV',0,1,20,'fatjet max subjet CSV','Events',999,-999,'fj1MaxCSV'))
plot.AddDistribution(root.Distribution('jet1CSV',0,1,20,'jet 1 CSV','Events',999,-999,'jet1CSV'))

recoil=None
if region=='signal':
  recoil=root.Distribution("pfmet",nRecoilBins,"MET [GeV]","a.u.")
elif 'single' in region: 
  recoil=root.Distribution('pfUWmag',nRecoilBins,'U(%s) [GeV]'%lep,"a.u.")
elif 'di' in region: 
  recoil=root.Distribution('pfUZmag',nRecoilBins,'U(%s%s) [GeV]'%(lep,lep),"a.u.")
elif region=='pho':
  recoil=root.Distribution('pfUAmag',nRecoilBins,'U(#gamma) [GeV]',"a.u.")
if recoil:
  setBins(recoil,recoilBins)
  plot.AddDistribution(recoil)
plot.AddDistribution(root.Distribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
plot.DrawAll(args.outdir+'/'+region+'_'+era+'_')
