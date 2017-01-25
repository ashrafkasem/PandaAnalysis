#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
args = parser.parse_args()
lumi = 36560.
blind=True
linear=False
sname = argv[0]

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.Monotop.MonojetSelection as sel
Load('Drawers','PlotUtility')

### DEFINE REGIONS ###

### LOAD PLOTTING UTILITY ###
plot = root.PlotUtility()
plot.DrawEmpty(True)
plot.Stack(True)
plot.FixRatio(3)
plot.Ratio(True)
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.SetDoOverflow()
plot.SetDoUnderflow()
plot.SetNormFactor(True)
plot.SetRatioLabel('#frac{SR-CR}{SR}')

plot.SetMCWeight('1')
plot.SetCut('1==1')

### DEFINE PROCESSES ###
sr     = root.Process('t#bar{t} (lep-veto)',root.kData); 
sr.additionalWeight = root.TCut(sel.weights['signal']%1)
sr.additionalCut = root.TCut(tAND(args.cut,sel.cuts['signal']))

cr     = root.Process('t#bar{t} (e or #mu)',root.kExtra1); 
cr.additionalWeight = root.TCut(sel.weights['singlemuon']%1)
cr.additionalCut = root.TCut(tAND(args.cut,tOR(sel.cuts['singlemuon'],sel.cuts['singleelectron'])))

processes = [sr,cr]

### ASSIGN FILES TO PROCESSES ###
sr.AddFile(baseDir+'TTbar.root')
cr.AddFile(baseDir+'TTbar.root')

for p in processes:
  plot.AddProcess(p)

plot.AddDistribution(root.Distribution('genWPlusPt',0,1000,20,'generated W^{+} p_{T} [GeV]','a.u./50 GeV'))
plot.AddDistribution(root.Distribution('genWMinusPt',0,1000,20,'generated W^{-} p_{T} [GeV]','a.u./50 GeV'))

plot.AddDistribution(root.Distribution('genWPlusEta', -5,5,20,'generated W^{+} #eta','a.u.'))
plot.AddDistribution(root.Distribution('genWMinusEta',-5,5,20,'generated W^{-} #eta','a.u.'))

plot.AddDistribution(root.Distribution('genTopPt',0,1000,20,'generated t p_{T} [GeV]','a.u./50 GeV'))
plot.AddDistribution(root.Distribution('genAntiTopPt',0,1000,20,'generated #bar{t} p_{T} [GeV]','a.u./50 GeV'))
plot.AddDistribution(root.Distribution('genTTPt',0,1000,20,'generated t#bar{t} p_{T} [GeV]','a.u./50 GeV'))

plot.AddDistribution(root.Distribution('genTopEta',-5,5,20,'generated t #eta','a.u.'))
plot.AddDistribution(root.Distribution('genAntiTopEta',-5,5,20,'generated #bar{t} #eta','a.u.'))
plot.AddDistribution(root.Distribution('genTTEta',-7,7,20,'generated t#bar{t} #eta','a.u.'))

plot.AddDistribution(root.Distribution('genTopEta-genAntiTopEta',-5,5,20,'#eta(gen t)-#eta(gen #bar{t})','a.u.'))
plot.AddDistribution(root.Distribution('TMath::Max(fabs(genTopEta),fabs(genAntiTopEta))',0,5,20,'max(#eta(gen t),#eta(gen #bar{t}))','a.u.',999,-999,'maxGenTopEta'))
plot.AddDistribution(root.Distribution('TMath::Max(fabs(genWPlusEta),fabs(genWMinusEta))',0,5,20,'max(#eta(gen W^{+}),#eta(gen W^{-}))','a.u.',999,-999,'maxGenWEta'))

plot.AddDistribution(root.Distribution('sf_tt',0.8,1.1,20,'13 TeV SF','a.u.'))
plot.AddDistribution(root.Distribution('sf_tt8TeV',0.5,1.3,20,'8 TeV SF','a.u.'))

plot.DrawAll(args.outdir+'/SRvCR_')
