#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
parser.add_argument('--region',metavar='region',type=str,default=None)
args = parser.parse_args()
lumi = 36560.
blind=True
region = args.region
sname = argv[0]

argv=[]
import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
#import PandaAnalysis.Monotop.NewPFSelection as sel
#import PandaAnalysis.Monotop.OldPFSelection as sel
#import PandaAnalysis.Monotop.NoTagPFSelection as sel
import PandaAnalysis.Monotop.NoMassPFSelection as sel
import PandaAnalysis.Tagging.cfg_v8 as cfg
Load('Drawers','PlotUtility')

### DEFINE REGIONS ###

cut = tAND(sel.cuts[args.region],args.cut)
PDebug(sname,'using cut: '+cut)

### LOAD PLOTTING UTILITY ###
plot = root.PlotUtility()
plot.Stack(True)
# plot.SetSignalScale(10)
plot.Ratio(True)
plot.FixRatio(0.6)
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.SetCut(root.TCut(cut))
plot.SetEvtNum("eventNumber")
if ('signal' in region) and blind:
  plot.SetEvtMod(5)
  plot.SetLumi(lumi/5000)
else:
  plot.SetLumi(lumi/1000)
plot.AddLumiLabel(True)
plot.SetDoOverflow()
plot.SetDoUnderflow()

weight = sel.weights[region]%lumi
plot.SetMCWeight(root.TCut(weight))
PDebug(sname,'using weight: '+weight)


### DEFINE PROCESSES ###
zjets     = root.Process('Z+jets',root.kZjets)
wjetsq     = root.Process('W+q',root.kWjets); wjetsq.additionalCut = root.TCut('abs(fj1HighestPtGen)!=21')
wjetsg     = root.Process('W+g',root.kExtra2); wjetsg.additionalCut = root.TCut('abs(fj1HighestPtGen)==21')
diboson   = root.Process('Diboson',root.kDiboson)
ttbar     = root.Process('Top [matched]',root.kTTbar); ttbar.additionalCut = root.TCut('(fj1IsMatched==1&&fj1GenSize<1.44)')
ttbarunmatched     = root.Process('Top [unmatched]',root.kExtra1); ttbarunmatched.additionalCut = root.TCut('(fj1IsMatched==0||fj1GenSize>1.44)')
#singletop = root.Process('Single t',root.kST)
qcd       = root.Process("QCD",root.kQCD)
gjetsq    = root.Process('#gamma+q',root.kGjets); gjetsq.additionalCut = root.TCut('abs(fj1HighestPtGen)!=21')
gjetsg    = root.Process('#gamma+g',root.kExtra3); gjetsg.additionalCut = root.TCut('abs(fj1HighestPtGen)==21')
data      = root.Process("Data",root.kData)
signal    = root.Process('m_{V}=1.7 TeV, m_{#chi}=100 GeV',root.kSignal)
#processes = [qcd,diboson,singletop,ttbar,wewk,zewk,wjets,zjets]
processes = [qcd,diboson,wjetsg,wjetsq,ttbarunmatched,ttbar,zjets]

### ASSIGN FILES TO PROCESSES ###
if 'signal' in region:
  zjets.AddFile(baseDir+'ZtoNuNu.root')
#  signal.AddFile(baseDir+'monotop-nr-v3-1700-100_med-1700_dm-100.root')
#  processes.append(signal)
else:
  zjets.AddFile(baseDir+'ZJets.root')
wjetsq.AddFile(baseDir+'WJets.root')
wjetsg.AddFile(baseDir+'WJets.root')
diboson.AddFile(baseDir+'Diboson.root')
ttbar.AddFile(baseDir+'TTbar.root')
ttbar.AddFile(baseDir+'SingleTop.root')
ttbarunmatched.AddFile(baseDir+'TTbar.root')
ttbarunmatched.AddFile(baseDir+'SingleTop.root')
if 'pho' in region:
  processes = [qcd,gjetsg,gjetsq]
  gjetsq.AddFile(baseDir+'GJets.root')
  gjetsg.AddFile(baseDir+'GJets.root')
  qcd.AddFile(baseDir+'SinglePhoton.root')
  qcd.additionalCut = root.TCut(sel.triggers['pho'])
  qcd.useCommonWeight = False
  qcd.additionalWeight = root.TCut('sf_phoPurity')
else:
  qcd.AddFile(baseDir+'QCD.root')

if any([x in region for x in ['singlemuonw','singleelectronw']]):
  processes = [qcd,diboson,zjets,ttbarunmatched,ttbar,wjetsg,wjetsq]
if any([x in region for x in ['singlemuontop','singleelectrontop']]):
  processes = [qcd,diboson,zjets,wjetsg,wjetsq,ttbarunmatched,ttbar]
if any([x in region for x in ['signal','muon']]):
  data.additionalCut = root.TCut(sel.triggers['met'])
  data.AddFile(baseDir+'MET.root')
  lep='#mu'
elif 'electron' in region:
  if 'di' in region:
    data.additionalCut = root.TCut(tOR(sel.triggers['ele'],sel.triggers['pho']))
  else:
    data.additionalCut = root.TCut(sel.triggers['ele'])
  data.AddFile(baseDir+'SingleElectron.root')
  lep='e'
elif region=='photon':
  data.additionalCut = root.TCut(sel.triggers['pho'])
  data.AddFile(baseDir+'SinglePhoton.root')
processes.append(data)

for p in processes:
  plot.AddProcess(p)

recoilBins = [250,280,310,350,400,450,600,1000]
nRecoilBins = len(recoilBins)-1

### CHOOSE DISTRIBUTIONS, LABELS ###

plot.AddDistribution(root.Distribution('fj1MSD',0,450,20,'fatjet m_{SD} [GeV]','Events/12.5 GeV'))

plot.AddDistribution(root.Distribution('top_ecf_bdt',-1,1,20,'ECF+#tau_{32}^{SD}+f_{rec} BDT','Events'))

plot.AddDistribution(root.Distribution('fj1Tau32SD',0,1,20,'Groomed #tau_{32}','Events',999,-999,'tau32SD'))

plot.AddDistribution(root.Distribution('fj1Tau32',0,1,20,'#tau_{32}','Events',999,-999,'tau32'))

plot.AddDistribution(root.Distribution('fj1HTTMass',40,450,20,'fatjet m_{HTT} [GeV]','Events/12.5 GeV'))

plot.AddDistribution(root.Distribution('fj1HTTFRec',0,1,20,'HTT f_{rec}','Events'))

plot.AddDistribution(root.Distribution('fj1ECFN_1_2_20/pow(fj1ECFN_1_2_10,2.00)',2,10,20,'e(1,2,2)/e(1,2,1)^{2}','Events',999,-999,'input0'))
plot.AddDistribution(root.Distribution('fj1ECFN_1_3_40/fj1ECFN_2_3_20',0,1,20,'e(1,3,4)/e(2,3,2)','Events',999,-999,'input1'))
plot.AddDistribution(root.Distribution('fj1ECFN_3_3_10/pow(fj1ECFN_1_3_40,.75)',.5,4,20,'e(3,3,1)/e(1,3,4)^{3/4}','Events',999,-999,'input2'))
plot.AddDistribution(root.Distribution('fj1ECFN_3_3_10/pow(fj1ECFN_2_3_20,.75)',0.4,1.4,20,'e(3,3,1)/e(2,3,2)^{3/4}','Events',999,-999,'input3'))
plot.AddDistribution(root.Distribution('fj1ECFN_3_3_20/pow(fj1ECFN_3_3_40,.5)',0,.25,20,'e(3,3,2)/e(3,3,4)^{1/2}','Events',999,-999,'input4'))
plot.AddDistribution(root.Distribution('fj1ECFN_1_4_20/pow(fj1ECFN_1_3_10,2)',0,2,20,'e(1,4,2)/e(1,3,1)^{2}','Events',999,-999,'input5'))
plot.AddDistribution(root.Distribution('fj1ECFN_1_4_40/pow(fj1ECFN_1_3_20,2)',0,2.5,20,'e(1,4,4)/e(1,3,2)^{2}','Events',999,-999,'input6'))
plot.AddDistribution(root.Distribution('fj1ECFN_2_4_05/pow(fj1ECFN_1_3_05,2)',1.25,2.5,20,'e(2,4,0.5)/e(1,3,0.5)^{2}','Events',999,-999,'input7'))
plot.AddDistribution(root.Distribution('fj1ECFN_2_4_10/pow(fj1ECFN_1_3_10,2)',1,4,20,'e(2,4,1)/e(1,3,1)^{2}','Events',999,-999,'input8'))
plot.AddDistribution(root.Distribution('fj1ECFN_2_4_10/pow(fj1ECFN_2_3_05,2)',0,1.5,20,'e(2,4,1)/e(2,3,0.5)^{2}','Events',999,-999,'input9'))
plot.AddDistribution(root.Distribution('fj1ECFN_2_4_20/pow(fj1ECFN_1_3_20,2)',0,5,20,'e(2,4,2)/e(1,3,2)^{2}','Events',999,-999,'input10'))

# if region=='signal':
plot.AddDistribution(root.Distribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
plot.DrawAll(args.outdir+'/'+region+'_')
