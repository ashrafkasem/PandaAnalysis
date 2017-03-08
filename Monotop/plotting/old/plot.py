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
lumi = 12900.
blind=False
linear=False
region = args.region
sname = argv[0]

argv=[]
import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
#import PandaAnalysis.Monotop.NoBTagSelection as sel
import PandaAnalysis.Monotop.NoTagSelection as sel
#import PandaAnalysis.Monotop.NewSelection as sel
#import PandaAnalysis.Monotop.OldSelection as sel
Load('Drawers','PlotUtility')

### DEFINE REGIONS ###

cut = tAND(sel.cuts[args.region],args.cut)
PInfo(sname,'using cut: '+cut)

### LOAD PLOTTING UTILITY ###
plot = root.PlotUtility()
plot.Stack(True)
plot.Logy(not(linear))
plot.SetLumi(lumi/1000)
# plot.SetSignalScale(10)
plot.Ratio(True)
plot.FixRatio(0.4)
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.AddLumiLabel(True)
plot.SetCut(root.TCut(cut))

weight = sel.weights[region]%lumi
plot.SetMCWeight(root.TCut(weight))
PInfo(sname,'using weight: '+weight)


### DEFINE PROCESSES ###
zjets     = root.Process('QCD Z+jets',root.kZjets)
wjets     = root.Process('QCD W+jets',root.kWjets)
diboson   = root.Process('Diboson',root.kDiboson)
ttbar     = root.Process('t#bar{t}',root.kTTbar)
# uttbar    = root.Process('t#bar{t} [unmatched]',root.kExtra1)
singletop = root.Process('Single t',root.kST)
qcd       = root.Process("QCD",root.kQCD)
gjets     = root.Process('#gamma+jets',root.kGjets)
data      = root.Process("Data",root.kData)
signal    = root.Process('m_{V}=1.7 TeV, m_{#chi}=100 GeV',root.kSignal)
#processes = [qcd,diboson,singletop,ttbar,wewk,zewk,wjets,zjets]
processes = [qcd,diboson,singletop,wjets,ttbar,zjets]

### ASSIGN FILES TO PROCESSES ###
if region=='signal':
  zjets.AddFile(baseDir+'ZtoNuNu.root')
  signal.AddFile(baseDir+'monotop-nr-v3-1700-100_med-1700_dm-100.root')
  processes.append(signal)
else:
  zjets.AddFile(baseDir+'ZJets.root')
wjets.AddFile(baseDir+'WJets.root')
diboson.AddFile(baseDir+'Diboson.root')
ttbar.AddFile(baseDir+'TTbar.root')
singletop.AddFile(baseDir+'SingleTop.root')
if 'pho' in region:
  processes = [qcd,gjets]
  gjets.AddFile(baseDir+'GJets.root')
  qcd.AddFile(baseDir+'SinglePhoton.root')
  qcd.additionalCut = root.TCut(sel.triggers['pho'])
  qcd.useCommonWeight = False
  qcd.additionalWeight = root.TCut('sf_phoPurity')
if any([x in region for x in ['signal','muon']]):
  data.additionalCut = root.TCut(sel.triggers['met'])
  PInfo(sname,'Using MET data')
  data.AddFile(baseDir+'MET.root')
  lep='#mu'
elif 'electron' in region:
  if 'di' in region:
    data.additionalCut = root.TCut(tOR(sel.triggers['pho'],sel.triggers['ele']))
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
if region in ['zmm','zee']:
  minval=0.1; maxval=5*10**3
else:
  minval=1; maxval=10**5
recoil=None
if region=='signal':
  recoil=root.Distribution("puppimet",nRecoilBins,"puppi MET [GeV]","Events/GeV")
  dphi=root.Distribution("dphipuppimet",0,3.14,20,"min#Delta#phi(jet,puppi MET)","Events")
elif any([x in region for x in ['singlemuonw','singleelectronw','singlemuontop','singleelectrontop']]):
  recoil=root.Distribution('UWmag',nRecoilBins,'U(%s) [GeV]'%lep,"Events/GeV")
  dphi=root.Distribution("dphiUW",0,3.14,20,"min#Delta#phi(jet,U(%s))"%lep,"Events")
  plot.AddDistribution(root.Distribution("dphipuppimet",0,3.14,20,"min#Delta#phi(jet,puppi MET)","Events"))
  plot.AddDistribution(root.Distribution("dphipfmet",0,3.14,20,"min#Delta#phi(jet,PF MET)","Events"))
  plot.AddDistribution(root.Distribution('looseLep1Pt',0,500,20,'Leading lep p_{T} [GeV]','Events/25 GeV'))
  plot.AddDistribution(root.Distribution('looseLep1Eta',-2.5,2.5,20,'Leading lep #eta','Events'))
elif any([x in region for x in ['dielectron','dimuon']]):
  recoil=root.Distribution('UZmag',nRecoilBins,'U(%s%s) [GeV]'%(lep,lep),"Events/GeV")
  dphi=root.Distribution("dphiUZ",0,3.14,20,"min#Delta#phi(jet,U(%s%s))"%(lep,lep),"Events")
  plot.AddDistribution(root.Distribution("dphipuppimet",0,3.14,20,"min#Delta#phi(jet,puppi MET)","Events"))
  plot.AddDistribution(root.Distribution("dphipfmet",0,3.14,20,"min#Delta#phi(jet,PF MET)","Events"))
  plot.AddDistribution(root.Distribution('looseLep1Pt',0,500,20,'Leading lep p_{T} [GeV]','Events/25 GeV'))
  plot.AddDistribution(root.Distribution('looseLep1Eta',-2.5,2.5,20,'Leading lep #eta','Events'))
  plot.AddDistribution(root.Distribution('looseLep2Pt',0,500,20,'Sub-leading lep p_{T} [GeV]','Events/25 GeV'))
  plot.AddDistribution(root.Distribution('looseLep2Eta',-2.5,2.5,20,'Sub-leading lep #eta','Events'))
elif region=='photon':
  recoil=root.Distribution('UAmag',nRecoilBins,'U(#gamma) [GeV]',"Events/GeV")
  dphi=root.Distribution("dphiUA",0,3.14,20,"min#Delta#phi(jet,U(#gamma))","Events")
  plot.AddDistribution(root.Distribution("dphipuppimet",0,3.14,20,"min#Delta#phi(jet,puppi MET)","Events"))
  plot.AddDistribution(root.Distribution("dphipfmet",0,3.14,20,"min#Delta#phi(jet,pf MET)","Events"))
  plot.AddDistribution(root.Distribution('loosePho1Pt',0,500,20,'Leading pho p_{T} [GeV]','Events/25 GeV'))
  plot.AddDistribution(root.Distribution('loosePho1Eta',-2.5,2.5,20,'Leading pho #eta','Events'))
  plot.AddDistribution(root.Distribution('loosePho1Phi',-3.142,3.142,20,'Leading pho #phi','Events'))
if recoil:
  setBins(recoil,recoilBins)
  plot.AddDistribution(recoil)
  plot.AddDistribution(dphi)

plot.AddDistribution(root.Distribution('top_ecf_bdt',-1,1,20,'ECF+#tau_{32}^{SD}+f_{rec} BDT','Events'))

plot.AddDistribution(root.Distribution('fj1Tau32SD',0,1,20,'Groomed #tau_{32}','Events',999,-999,'tau32SD'))

plot.AddDistribution(root.Distribution('fj1Tau32',0,1,20,'#tau_{32}','Events',999,-999,'tau32'))

plot.AddDistribution(root.Distribution('jet1Pt',15,500,20,'leading jet p_{T} [GeV]','Events'))

# plot.AddDistribution(root.Distribution('nJet',-0.5,8.5,9,'N_{jet}','Events'))

plot.AddDistribution(root.Distribution('fj1Pt',250,1000,20,'fatjet p_{T} [GeV]','Events/37.5 GeV'))

plot.AddDistribution(root.Distribution('fj1MSD',40,450,20,'fatjet m_{SD} [GeV]','Events/12.5 GeV'))

plot.AddDistribution(root.Distribution('max(0.01,fj1MaxCSV)',0,1,20,'fatjet max subjet CSV','Events',999,-999,'fj1MaxCSV'))

plot.AddDistribution(root.Distribution('max(0.01,jet1CSV)',0,1,20,'jet 1 CSV','Events',999,-999,'jet1CSV'))

'''
'''

# plot.AddDistribution(root.Distribution('fj1HTTMass',40,450,20,'fatjet m_{HTT} [GeV]','Events/12.5 GeV'))

# plot.AddDistribution(root.Distribution('fj1HTTFRec',0,1,20,'HTT f_{rec}','Events'))

# plot.AddDistribution(root.Distribution('fj1ECFN_1_2_20/pow(fj1ECFN_1_2_10,2.00)',2,10,20,'e(1,2,2)/e(1,2,1)^{2}','Events',999,-999,'input0'))
# plot.AddDistribution(root.Distribution('fj1ECFN_1_3_40/fj1ECFN_2_3_20',0,1,20,'e(1,3,4)/e(2,3,2)','Events',999,-999,'input1'))
# plot.AddDistribution(root.Distribution('fj1ECFN_3_3_10/pow(fj1ECFN_1_3_40,.75)',.5,4,20,'e(3,3,1)/e(1,3,4)^{3/4}','Events',999,-999,'input2'))
# plot.AddDistribution(root.Distribution('fj1ECFN_3_3_10/pow(fj1ECFN_2_3_20,.75)',0.4,1.4,20,'e(3,3,1)/e(2,3,2)^{3/4}','Events',999,-999,'input3'))
# plot.AddDistribution(root.Distribution('fj1ECFN_3_3_20/pow(fj1ECFN_3_3_40,.5)',0,.25,20,'e(3,3,2)/e(3,3,4)^{1/2}','Events',999,-999,'input4'))
# plot.AddDistribution(root.Distribution('fj1ECFN_1_4_20/pow(fj1ECFN_1_3_10,2)',0,2,20,'e(1,4,2)/e(1,3,1)^{2}','Events',999,-999,'input5'))
# plot.AddDistribution(root.Distribution('fj1ECFN_1_4_40/pow(fj1ECFN_1_3_20,2)',0,2.5,20,'e(1,4,4)/e(1,3,2)^{2}','Events',999,-999,'input6'))
# plot.AddDistribution(root.Distribution('fj1ECFN_2_4_05/pow(fj1ECFN_1_3_05,2)',1.25,2.5,20,'e(2,4,0.5)/e(1,3,0.5)^{2}','Events',999,-999,'input7'))
# plot.AddDistribution(root.Distribution('fj1ECFN_2_4_10/pow(fj1ECFN_1_3_10,2)',1,4,20,'e(2,4,1)/e(1,3,1)^{2}','Events',999,-999,'input8'))
# plot.AddDistribution(root.Distribution('fj1ECFN_2_4_10/pow(fj1ECFN_2_3_05,2)',0,1.5,20,'e(2,4,1)/e(2,3,0.5)^{2}','Events',999,-999,'input9'))
# plot.AddDistribution(root.Distribution('fj1ECFN_2_4_20/pow(fj1ECFN_1_3_20,2)',0,5,20,'e(2,4,2)/e(1,3,2)^{2}','Events',999,-999,'input10'))

# if region=='signal':
plot.AddDistribution(root.Distribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
plot.DrawAll(args.outdir+'/'+region+'_')
