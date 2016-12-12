#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--region',metavar='region',type=str,default=None)
parser.add_argument('--ptbin',metavar='ptbin',type=int,default=0)
args = parser.parse_args()
lumi = 36560.
blind=True
linear=False
region = args.region
sname = argv[0]

ptbins = [
    (0,1000),
    (0,50),
    (50,100),
    (100,150),
    (150,200),
    (200,1000),
    ]

argv=[]
import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.Monotop.CSVSelection as sel
Load('Drawers','PlotUtility')

### DEFINE REGIONS ###

ptcut = 'isojet1Pt>%f && isojet1Pt<%f'%(ptbins[args.ptbin])
cut = tAND(ptcut,sel.cuts[args.region])
PInfo(sname,'using cut: '+cut)

### LOAD PLOTTING UTILITY ###
plot = root.PlotUtility()
plot.Stack(True)
plot.Logy(False)
plot.Ratio(True)
plot.FixRatio(0.4)
plot.SetTDRStyle()
plot.InitLegend()
plot.SetLumi(lumi/1000.)
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.SetCut(root.TCut(cut))
plot.AddLumiLabel(True)
plot.SetDoOverflow()
plot.SetDoUnderflow()
plot.SetNormFactor(True)

weight = sel.weights[region]%lumi
plot.SetMCWeight(root.TCut(weight))
PInfo(sname,'using weight: '+weight)

### DEFINE PROCESSES ###
zjets     = root.Process('Z+jets',root.kZjets)
wjets     = root.Process('W+jets',root.kWjets)
diboson   = root.Process('Diboson',root.kDiboson)
ttbar     = root.Process('Top (HF)',root.kTTbar); ttbar.additionalCut = root.TCut('isojet1Flav!=0')
ttbarlf   = root.Process('Top (LF)',root.kExtra1); ttbarlf.additionalCut = root.TCut('isojet1Flav==0')
qcd       = root.Process("QCD",root.kQCD)
gjets     = root.Process('#gamma+jets (LF)',root.kGjets); gjets.additionalCut = root.TCut('isojet1Flav==0')
gjetshf   = root.Process('#gamma+jets (HF)',root.kExtra4); gjetshf.additionalCut = root.TCut('isojet1Flav!=0')
data      = root.Process("Data",root.kData)
signal    = root.Process('m_{V}=1.7 TeV, m_{#chi}=100 GeV',root.kSignal)
processes = [qcd,diboson,zjets,wjets,ttbarlf,ttbar]

### ASSIGN FILES TO PROCESSES ###
zjets.AddFile(baseDir+'ZJets.root')
wjets.AddFile(baseDir+'WJets.root')
diboson.AddFile(baseDir+'Diboson.root')
ttbar.AddFile(baseDir+'TTbar.root')
ttbarlf.AddFile(baseDir+'TTbar.root')
ttbar.AddFile(baseDir+'SingleTop.root')
ttbarlf.AddFile(baseDir+'SingleTop.root')
if 'pho' in region:
  processes = [qcd,gjetshf,gjets]
  gjets.AddFile(baseDir+'GJets.root')
  gjetshf.AddFile(baseDir+'GJets.root')
  qcd.AddFile(baseDir+'SinglePhoton.root')
  qcd.additionalCut = root.TCut(sel.triggers['pho'])
  qcd.useCommonWeight = False
  qcd.additionalWeight = root.TCut('sf_phoPurity')
  data.additionalCut = root.TCut(sel.triggers['pho'])
  data.AddFile(baseDir+'SinglePhoton.root')
else:
  qcd.AddFile(baseDir+'QCD.root')
  data.additionalCut = root.TCut(tOR(sel.triggers['ele'],sel.triggers['met']))
  data.AddFile(baseDir+'METSingleElectron.root')

processes.append(data)

for p in processes:
  plot.AddProcess(p)

recoilBins = [250,280,310,350,400,450,600,1000]
nRecoilBins = len(recoilBins)-1

### CHOOSE DISTRIBUTIONS, LABELS ###
recoil=None
lep = 'l'
if  'ttbar' in region:
  recoil=root.Distribution('pfUWmag',nRecoilBins,'U(%s) [GeV]'%lep,"a.u./GeV")
  dphi=root.Distribution("dphiUW",0,3.14,20,"min#Delta#phi(jet,U(%s))"%lep,"a.u.")
  plot.AddDistribution(root.Distribution("dphipfmet",0,3.14,20,"min#Delta#phi(jet,PF MET)","a.u."))
  plot.AddDistribution(root.Distribution('looseLep1Pt',0,500,20,'Leading lep p_{T} [GeV]','a.u./25 GeV'))
  plot.AddDistribution(root.Distribution('looseLep1Eta',-2.5,2.5,20,'Leading lep #eta','a.u.'))
  plot.AddDistribution(root.Distribution('looseLep2Pt',0,500,20,'Sub-leading lep p_{T} [GeV]','a.u./25 GeV'))
  plot.AddDistribution(root.Distribution('looseLep2Eta',-2.5,2.5,20,'Sub-leading lep #eta','a.u.'))
elif 'photon' in region:
  recoil=root.Distribution('pfUAmag',nRecoilBins,'U(#gamma) [GeV]',"a.u./GeV")
  dphi=root.Distribution("dphiUA",0,3.14,20,"min#Delta#phi(jet,U(#gamma))","a.u.")
  plot.AddDistribution(root.Distribution("dphipfmet",0,3.14,20,"min#Delta#phi(jet,pf MET)","a.u."))
  plot.AddDistribution(root.Distribution('loosePho1Pt',175,1175,20,'Leading pho p_{T} [GeV]','a.u./50 GeV'))
  plot.AddDistribution(root.Distribution('loosePho1Eta',-2.5,2.5,20,'Leading pho #eta','a.u.'))
  plot.AddDistribution(root.Distribution('loosePho1Phi',-3.142,3.142,20,'Leading pho #phi','a.u.'))

if recoil:
  setBins(recoil,recoilBins)
  plot.AddDistribution(recoil)
  plot.AddDistribution(dphi)


plot.AddDistribution(root.Distribution('jet1Pt',15,1000,20,'leading jet p_{T} [GeV]','a.u.'))

plot.AddDistribution(root.Distribution('isojet1Pt',15,1000,20,'leading isojet p_{T} [GeV]','a.u.'))

plot.AddDistribution(root.Distribution('nJet',-0.5,8.5,9,'N_{jet}','a.u.'))

plot.AddDistribution(root.Distribution('jet1CSV',0,1,10,'jet 1 CSV','a.u.',999,-999,'jet1CSV'))

plot.AddDistribution(root.Distribution('isojet1CSV',0,1,10,'isojet 1 CSV','a.u.',999,-999,'isojet1CSV'))

plot.AddDistribution(root.Distribution('jet2CSV',0,1,10,'jet 2 CSV','a.u.',999,-999,'jet2CSV'))

plot.AddDistribution(root.Distribution('isojet2CSV',0,1,10,'isojet 2 CSV','a.u.',999,-999,'isojet2CSV'))

# if region=='signal':
plot.AddDistribution(root.Distribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
if args.ptbin==0:
  plot.DrawAll(args.outdir+'/'+region+'_')
else:
  plot.DrawAll(args.outdir+'/pt%i/'%(args.ptbin)+region+'_')
