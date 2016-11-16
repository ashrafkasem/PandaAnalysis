#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_ZEYNEPDIR')+'/merged/' 
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--which',metavar='which',type=str)
parser.add_argument('--region',metavar='region',type=str,default="signal")
args = parser.parse_args()
lumi = 12900.
blind=True
linear=False
region = args.region
sname = argv[0]

argv=[]
import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.VBF.LooseSelection as sel
Load('Drawers','PlotUtility')

### DEFINE REGIONS ###

cut = sel.cuts[args.region]
if args.which=='jot1Eta':
  cut=removeCut(removeCut(cut,'leadingJet_outaccp'),'jet1isMonoJetIdNew')
  label=''
elif args.which=='dphi_only':
  cut=removeCut(cut,'fabs(minJetMetDPhi_withendcap)')
  label='only_'
else:
  label=''
  cut=removeCut(removeCut(removeCut(cut,'leadingJet_outaccp'),'jet1isMonoJetIdNew'),'fabs(minJetMetDPhi_withendcap)')
PInfo(sname,'using cut: '+cut)

### LOAD PLOTTING UTILITY ###
plot = root.PlotUtility()
plot.Stack(True)
plot.Logy(not(linear))
if 'signal' in region and blind:
  plot.SetLumi(lumi/5000.)
else:
  plot.SetLumi(lumi/1000)
plot.SetSignalScale(10)
plot.Ratio(True)
# plot.FixRatio()
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.AddLumiLabel(True)
plot.SetEvtNum("eventNum")
if region=='signal' and blind:
  plot.SetEvtMod(5)
plot.SetCut(root.TCut(cut))

weight = '%f*%s'%(lumi,sel.weights[region])
plot.SetMCWeight(root.TCut(weight))
PInfo(sname,'using weight: '+weight)


### DEFINE PROCESSES ###
zjets     = root.Process('QCD Z+jets',root.kZjets); zjets.additionalWeight = root.TCut('zkfactor*ewk_z')
zewk      = root.Process('EWK Z+jets',root.kExtra2)
wjets     = root.Process('QCD W+jets',root.kWjets); wjets.additionalWeight = root.TCut('ewk_w')
wewk      = root.Process('EWK W+jets',root.kExtra3)
diboson   = root.Process('Diboson',root.kDiboson)
ttbar     = root.Process('t#bar{t}',root.kTTbar)
singletop = root.Process('Single t',root.kST)
qcd       = root.Process("QCD",root.kQCD)
gjets     = root.Process('#gamma+jets',root.kGjets); gjets.additionalWeight = root.TCut('akfactor*ewk_a')
vbf       = root.Process("VBF H#rightarrowInv",root.kSignal)
data      = root.Process("Data",root.kData)
#processes = [qcd,diboson,singletop,ttbar,wewk,zewk,wjets,zjets]
processes = [diboson,singletop,ttbar,wewk,zewk,wjets,zjets]

### ASSIGN FILES TO PROCESSES ###
if region=='signal':
  zjets.AddFile(baseDir+'ZtoNuNu.root')
  zewk.AddFile(baseDir+'EWKZtoNuNu.root')
else:
  zjets.AddFile(baseDir+'ZJets.root')
  zewk.AddFile(baseDir+'EWKZJets.root')
wjets.AddFile(baseDir+'WJets.root')
wewk.AddFile(baseDir+'EWKWJets.root')
diboson.AddFile(baseDir+'Diboson.root')
ttbar.AddFile(baseDir+'TTbar.root')
singletop.AddFile(baseDir+'SingleTop.root')
qcd.AddFile(baseDir+'QCD.root')
if 'pho' in region:
  processes = [qcd,gjets]
  gjets.AddFile(baseDir+'GJets.root')

if 'signal' in region:
  vbf.AddFile(baseDir+'VBF_H125.root')
  processes += [vbf]
if region in  ['signal','zmm','wmn']:
  data.additionalCut = root.TCut(tAND(sel.triggers['met'],'runNum<=276811'))
  data.AddFile(baseDir+'MET.root')
  lep='#mu'
elif region in ['zee','wen']:
  data.additionalCut = root.TCut(sel.triggers['ele'])
  data.AddFile(baseDir+'SingleElectron.root')
  lep='e'
elif region=='pho':
  data.additionalCut = root.TCut(sel.triggers['pho'])
  data.AddFile(baseDir+'SinglePhoton.root')
processes.append(data)

for p in processes:
  plot.AddProcess(p)

#recoilBins = [200,250,300,350,400,500,600,750,1000]
recoilBins = [200,250,300,350,400,500,600,1000]
nRecoilBins = len(recoilBins)-1

### CHOOSE DISTRIBUTIONS, LABELS ###
# plot.AddDistribution(root.Distribution("1",0,2,1,"dummy","dummy"))
if region in ['zmm','zee']:
  minval=0.1; maxval=5*10**3
else:
  minval=1; maxval=10**5
if args.which=='jjDEta':
  plot.AddDistribution(root.Distribution("jot1Eta",-5,5,20,"Jet 1 #eta","Events",minval,maxval))
else:
  plot.AddDistribution(root.Distribution('fabs(minJetMetDPhi_withendcap)',0,3.142,20,'min #Delta #phi(jet, MET)','Events',minval,maxval,'minDPhiJetMet'))

### DRAW AND CATALOGUE ###
plot.DrawAll(args.outdir+'/'+label+region+'_')
