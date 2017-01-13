#!/usr/bin/env python

from os import system,getenv
from sys import argv,exit
import argparse

### SET GLOBAL VARIABLES ###
rerecoDir = getenv('PANDA_ZEYNEPDIR')+'/merged/' 
promptDir = getenv('PANDA_ZEYNEPDIR_PROMPT')+'/merged/' 
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
#import PandaAnalysis.VBF.MonojetSelection as sel
import PandaAnalysis.VBF.LooseSelection as sel
Load('Drawers','PlotUtility')

all_runs = {
    'B' : (272007,275376),
    'C' : (275657,276283),
    'D' : (276315,276811),
    'E' : (276831,277420),
    'F' : (277772,278808),
    'G' : (278820,280385),
    'H' : (280919,284044),
    }
runs = all_runs[era]
cut = tAND(tAND(sel.cuts[args.region],args.cut),'runNum>%i && runNum<%i'%(runs[0],runs[1]))

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
plot.SetEvtNum("eventNum")
if region=='signal':
  cut = tAND(cut,'eventNum%5==0')
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
if region in  ['signal','zmm','wmn']:
  rereco.AddFile(rerecoDir+'MET.root')
  prompt.AddFile(promptDir+'MET.root')
  lep='#mu'
elif region in ['zee','wen']:
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
if not (region in ['signal','pho']):
  plot.AddDistribution(root.Distribution("lep1Pt",0,500,20,"lep 1 p_{T} [GeV]","a.u."))
  plot.AddDistribution(root.Distribution("lep1Eta",-2.5,2.5,20,"lep 1 #eta","Events"))
if region in ['pho']:
  plot.AddDistribution(root.Distribution("photonPt",175,675,20,"photon 1 p_{T} [GeV]","a.u."))
  plot.AddDistribution(root.Distribution("photonEta",-1.4442,1.4442,20,"photon 1 #eta","Events"))
if region in ['wmn','wen']:
  plot.AddDistribution(root.Distribution('trueMet',0,500,20,'true MET [GeV]','a.u.'))
  plot.AddDistribution(root.Distribution('mt',0,160,20,'M_{T} [GeV]','a.u.'))
plot.AddDistribution(root.Distribution('fabs(SignedDeltaPhi(jot1Phi,metPhi))',0,3.142,20,'#Delta #phi(jet 1, MET)','Events',999,-999,'jet1DPhiMet'))
plot.AddDistribution(root.Distribution('fabs(SignedDeltaPhi(jot2Phi,metPhi))',0,3.142,20,'#Delta #phi(jet 1, MET)','Events',999,-999,'jet2DPhiMet'))
plot.AddDistribution(root.Distribution('fabs(minJetMetDPhi_withendcap)',0,3.142,20,'min #Delta #phi(jet, MET)','Events',999,-999,'minDPhiJetMet'))
plot.AddDistribution(root.Distribution('fixed_mjj',0,4000,20,"Dijet mass [GeV]","a.u."))
plot.AddDistribution(root.Distribution("jjDEta",0,10,20,"Delta #eta leading jets","Events"))
plot.AddDistribution(root.Distribution("jot1Pt",40,760,24,"Jet 1 p_{T} [GeV]","a.u."))
plot.AddDistribution(root.Distribution("jot2Pt",40,760,24,"Jet 2 p_{T} [GeV]","a.u."))
plot.AddDistribution(root.Distribution("jot1Eta",-5,5,20,"Jet 1 #eta","Events"))
plot.AddDistribution(root.Distribution("jot2Eta",-5,5,20,"Jet 2 #eta","Events"))
plot.AddDistribution(root.Distribution("fabs(SignedDeltaPhi(jot1Phi,jot2Phi))",0,3.142,20,"#Delta #phi leading jets","Events",999,-999,'jjDPhi'))

recoil=None
if region=='signal':
  recoil=root.Distribution("met",nRecoilBins,"MET [GeV]","a.u.")
elif any([x in region for x in ['wen','wmn']]):
  recoil=root.Distribution('met',nRecoilBins,'U(%s) [GeV]'%lep,"a.u.")
elif any([x in region for x in ['zee','zmm']]):
  recoil=root.Distribution('met',nRecoilBins,'U(%s%s) [GeV]'%(lep,lep),"a.u.")
elif region=='pho':
  recoil=root.Distribution('met',nRecoilBins,'U(#gamma) [GeV]',"a.u.")
if recoil:
  setBins(recoil,recoilBins)
  plot.AddDistribution(recoil)
plot.AddDistribution(root.Distribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
plot.DrawAll(args.outdir+'/'+region+'_'+era+'_')
