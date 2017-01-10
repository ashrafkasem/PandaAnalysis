#!/usr/bin/env python

from os import system,getenv
from sys import argv,exit
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_ZEYNEPDIR')+'/merged/' 
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
parser.add_argument('--region',metavar='region',type=str,default=None)
args = parser.parse_args()
lumi = 35000
blind=True
linear=False
region = args.region
sname = argv[0]

argv=[]
import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
#import PandaAnalysis.VBF.TightSelection as sel
#import PandaAnalysis.VBF.Selection as sel
import PandaAnalysis.VBF.LooseSelection as sel
#import PandaAnalysis.VBF.MonojetSelection as sel
Load('Drawers','PlotUtility')

### DEFINE REGIONS ###

cut = tAND(sel.cuts[args.region],args.cut)
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
plot.FixRatio(0.4)
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
wjets     = root.Process('QCD W+jets',root.kWjets); wjets.additionalWeight = root.TCut('wkfactor*ewk_w')
wewk      = root.Process('EWK W+jets',root.kExtra3)
diboson   = root.Process('Diboson',root.kDiboson)
ttbar     = root.Process('t#bar{t}',root.kTTbar)
singletop = root.Process('Single t',root.kST)
qcd       = root.Process("QCD",root.kQCD)
gjets     = root.Process('#gamma+jets',root.kGjets); gjets.additionalWeight = root.TCut('akfactor*ewk_a')
vbf       = root.Process("VBF H#rightarrowInv",root.kSignal)
data      = root.Process("Data",root.kData)
processes = [qcd,diboson,singletop,ttbar,wewk,zewk,wjets,zjets]
#processes = [diboson,singletop,ttbar,wewk,zewk,wjets,zjets]

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
  data.additionalCut = root.TCut(sel.triggers['met'])
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
if region in ['zmm','zee']:
  minval=0.1; maxval=5*10**3
else:
  minval=1; maxval=10**5
if not (region in ['signal','pho']):
  plot.AddDistribution(root.Distribution("lep1Pt",0,500,20,"lep 1 p_{T} [GeV]","Events/25 GeV",minval,maxval))
  plot.AddDistribution(root.Distribution("lep1Eta",-2.5,2.5,20,"lep 1 #eta","Events",minval,maxval))
if region in ['pho']:
  plot.AddDistribution(root.Distribution("photonPt",175,675,20,"photon 1 p_{T} [GeV]","Events/25 GeV",minval,maxval))
  plot.AddDistribution(root.Distribution("photonEta",-1.4442,1.4442,20,"photon 1 #eta","Events",minval,maxval))
if region in ['wmn','wen']:
  plot.AddDistribution(root.Distribution('trueMet',0,500,20,'true MET [GeV]','Events/25 GeV',minval,maxval))
  plot.AddDistribution(root.Distribution('mt',0,160,20,'M_{T} [GeV]','Events/8 GeV',minval,maxval))
plot.AddDistribution(root.Distribution('fabs(jet1DPhiMet)',0,3.142,20,'#Delta #phi(jet 1, MET)','Events',minval,maxval,'jet1DPhiMet'))
plot.AddDistribution(root.Distribution('fabs(minJetMetDPhi_withendcap)',0,3.142,20,'min #Delta #phi(jet, MET)','Events',minval,maxval,'minDPhiJetMet'))
plot.AddDistribution(root.Distribution(sel.fixedmjj,0,4000,20,"Dijet mass [GeV]","Events/200 GeV",minval,maxval))
plot.AddDistribution(root.Distribution("jjDEta",0,10,20,"Delta #eta leading jets","Events",minval,maxval))
plot.AddDistribution(root.Distribution("jot1Pt",40,760,24,"Jet 1 p_{T} [GeV]","Events/30 GeV",minval,maxval))
plot.AddDistribution(root.Distribution("jot2Pt",40,760,24,"Jet 2 p_{T} [GeV]","Events/30 GeV",minval,maxval))
plot.AddDistribution(root.Distribution("jot1Eta",-5,5,20,"Jet 2 #eta","Events",minval,maxval))
plot.AddDistribution(root.Distribution("jot2Eta",-5,5,20,"Jet 2 #eta","Events",minval,maxval))
plot.AddDistribution(root.Distribution("fabs(SignedDeltaPhi(jot1Phi,jot2Phi))",0,3.142,20,"#Delta #phi leading jets","Events",minval,maxval,'jjDPhi'))

recoil=None
if region=='signal':
  recoil=root.Distribution("met",nRecoilBins,"MET [GeV]","Events/GeV")
elif any([x in region for x in ['wen','wmn']]):
  recoil=root.Distribution('met',nRecoilBins,'U(%s) [GeV]'%lep,"Events/GeV")
elif any([x in region for x in ['zee','zmm']]):
  recoil=root.Distribution('met',nRecoilBins,'U(%s%s) [GeV]'%(lep,lep),"Events/GeV")
elif region=='pho':
  recoil=root.Distribution('met',nRecoilBins,'U(#gamma) [GeV]',"Events/GeV")
if recoil:
  setBins(recoil,recoilBins)
  plot.AddDistribution(recoil)
plot.AddDistribution(root.Distribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
plot.DrawAll(args.outdir+'/'+region+'_')
