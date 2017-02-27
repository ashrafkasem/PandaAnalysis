#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--basedir',metavar='basedir',type=str,default=None)
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
parser.add_argument('--region',metavar='region',type=str,default=None)
parser.add_argument('--tt',metavar='tt',type=str,default='')
args = parser.parse_args()
lumi = 36560.
blind=True
linear=False
region = args.region
sname = argv[0]
if args.basedir:
	baseDir = args.basedir

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.Monotop.TestSelection as sel
from PandaCore.Drawers.plot_utility import *

### DEFINE REGIONS ###

cut = tAND(sel.cuts[args.region],args.cut)

PInfo(sname,'using cut: '+cut)

### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.Stack(True)
plot.Logy(not(linear))
# plot.SetSignalScale(10)
plot.Ratio(True)
plot.FixRatio(0.4)
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.cut = cut
plot.SetEvtNum("eventNumber")
if ('signal' in region) and blind:
  plot.SetEvtMod(5)
  plot.SetLumi(lumi/5000)
  plot.AddPlotLabel('Every 5th event',.18,.7,False,42,.04)
else:
  plot.SetLumi(lumi/1000)
plot.AddLumiLabel(True)
plot.do_overflow = True
plot.do_underflow = True

weight = sel.weights[region]%lumi
plot.mc_weight = weight

plot.add_systematic('QCD scale','scaleUp','scaleDown',root.kRed+2)
plot.add_systematic('PDF','pdfUp','pdfDown',root.kBlue+2)

### DEFINE PROCESSES ###
zjets     = Process('Z+jets',root.kZjets)
wjets     = Process('W+jets',root.kWjets)
diboson   = Process('Diboson',root.kDiboson)
ttbar     = Process('t#bar{t}',root.kTTbar); ttbar.additional_weight = '730/831.'
singletop = Process('Single t',root.kST)
qcd       = Process("QCD",root.kQCD)
gjets     = Process('#gamma+jets',root.kGjets)
data      = Process("Data",root.kData)
signal    = Process('m_{V}=1.7 TeV, m_{#chi}=100 GeV',root.kSignal)
#processes = [qcd,diboson,singletop,ttbar,wewk,zewk,wjets,zjets]
processes = [qcd,diboson,singletop,wjets,ttbar,zjets]

### ASSIGN FILES TO PROCESSES ###
if 'signal' in region:
  zjets.add_file(baseDir+'ZtoNuNu.root')
else:
  zjets.add_file(baseDir+'ZJets.root')
wjets.add_file(baseDir+'WJets.root')
diboson.add_file(baseDir+'Diboson.root')
ttbar.add_file(baseDir+'TTbar.root')
singletop.add_file(baseDir+'SingleTop.root')
if 'pho' in region:
  processes = [qcd,gjets]
  gjets.add_file(baseDir+'GJets.root')
  qcd.add_file(baseDir+'SinglePhoton.root')
  qcd.additional_cut = sel.triggers['pho']
  qcd.use_common_weight = False
  qcd.additional_weight = 'sf_phoPurity'
else:
  qcd.add_file(baseDir+'QCD.root')

if any([x in region for x in ['singlemuonw','singleelectronw']]):
  processes = [qcd,diboson,singletop,zjets,ttbar,wjets,]
if any([x in region for x in ['singlemuontop','singleelectrontop']]):
  processes = [qcd,diboson,singletop,zjets,wjets,ttbar]
if any([x in region for x in ['signal','muon']]):
  data.additional_cut = sel.triggers['met']
  PInfo(sname,'Using MET data')
  data.add_file(baseDir+'MET.root')
  lep='#mu'
elif 'electron' in region:
  if 'di' in region:
    data.additional_cut = tOR(sel.triggers['ele'],sel.triggers['pho'])
  else:
    data.additional_cut = sel.triggers['ele']
  data.add_file(baseDir+'SingleElectron.root')
  lep='e'
elif region=='photon':
  data.additional_cut = sel.triggers['pho']
  data.add_file(baseDir+'SinglePhoton.root')
processes.append(data)

for p in processes:
  plot.add_process(p)

recoilBins = [250,280,310,350,400,450,600,1000]
nRecoilBins = len(recoilBins)-1

### CHOOSE DISTRIBUTIONS, LABELS ###
if 'signal' in region:
  recoil=VDistribution("pfmet",recoilBins,"PF MET [GeV]","Events/GeV")
elif any([x in region for x in ['singlemuonw','singleelectronw','singlemuontop','singleelectrontop','singlemuon','singleelectron']]):
  recoil=VDistribution("pfUWmag",recoilBins,"PF U(%s) [GeV]"%(lep),"Events/GeV")
elif any([x in region for x in ['dielectron','dimuon']]):
  recoil=VDistribution("pfUZmag",recoilBins,"PF U(%s%s) [GeV]"%(lep,lep),"Events/GeV")
elif region=='photon':
  recoil=VDistribution("pfUAmag",recoilBins,"PF U(#gamma) [GeV]","Events/GeV")

plot.add_distribution(recoil)

plot.add_distribution(FDistribution('fj1MSD',0,400,20,'fatjet m_{SD} [GeV]','Events'))
plot.add_distribution(FDistribution('fj1Pt',250,600,20,'fatjet p_{T} [GeV]','Events'))
plot.add_distribution(FDistribution('top_ecf_bdt',0,1,20,'Top BDT','Events'))
plot.add_distribution(FDistribution('fj1MaxCSV',0,1,20,'fatjet max CSV','Events'))
plot.add_distribution(FDistribution('jet1CSV',0,1,20,'jet 1 CSV','Events',filename='jet1CSV'))
plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
plot.draw_all(args.outdir+'/'+region+'_')
