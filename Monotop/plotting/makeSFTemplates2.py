#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv,exit
import argparse

basedir = getenv('PANDA_FLATDIR')
figsdir = basedir+'/figs'

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--indir',metavar='indir',type=str,default=basedir)
parser.add_argument('--outdir',metavar='outdir',type=str,default=figsdir)
parser.add_argument('--disc',metavar='disc',type=str,default='top_ecfv8_bdt')
parser.add_argument('--tagged',metavar='tagged',type=str,default='True')
parser.add_argument('--cutlo',metavar='cutlo',type=float,default=0.1)
parser.add_argument('--cuthi',metavar='cuthi',type=float,default=0.45)
parser.add_argument('--pt',metavar='pt',type=str,default='inc')
args = parser.parse_args()

figsdir = args.outdir+'/'+args.pt
basedir = args.indir
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
import PandaAnalysis.Monotop.NoMassPFSelection as sel
tcut = root.TCut
Load('Drawers','PlotUtility')

### SET GLOBAL VARIABLES ###
lumi = 36560.
logy=False

cut = tOR(sel.cuts['singlemuontop'],sel.cuts['singleelectrontop'])
weight = sel.weights['singlemuontop']%lumi
label = ''

disclabels = {
    'fj1Tau32SD':'#tau_{32}^{SD}',
    'top_ecfv6_bdt':'ECF BDT',
    'top_ecfv7_bdt':'ECF+#tau_{32}^{SD} BDT',
    'top_ecfv8_bdt':'ECF+#tau_{32}^{SD}+f_{rec} BDT',
    'top_ecf_bdt':'ECF+#tau_{32}^{SD}+f_{rec} BDT',
    }

label += args.disc+'_'
if args.tagged.upper()=='TRUE':
  label += 'pass_'
  cut = tAND(cut,'%f<%s && %s<%f'%(args.cutlo,args.disc,args.disc,args.cuthi))
  plotlabel = '%.2f<%s<%.2f'%(args.cutlo,disclabels[args.disc],args.cuthi)
else:
  label += 'fail_'
  cut = tAND(cut,'!(%f<%s && %s<%f)'%(args.cutlo,args.disc,args.disc,args.cuthi))
  plotlabel = '!(%.2f<%s<%.2f)'%(args.cutlo,disclabels[args.disc],args.cuthi)

if args.pt=='' or args.pt=='inc':
  cut = tAND(cut,'fj1Pt>250 && fj1Pt<1000')
  plotlabel = '#splitline{%s}{250 < p_{T} < 1000 GeV}'%plotlabel
elif args.pt=='lo':
  cut = tAND(cut,'fj1Pt>250 && fj1Pt<450')
  plotlabel = '#splitline{%s}{250 < p_{T} < 450 GeV}'%plotlabel
elif args.pt=='med':
  cut = tAND(cut,'fj1Pt>450 && fj1Pt<1000')
  plotlabel = '#splitline{%s}{450 < p_{T} < 500 GeV}'%plotlabel
elif args.pt=='hi':
  cut = tAND(cut,'fj1Pt>500 && fj1Pt<1000')
  plotlabel = '#splitline{%s}{500 < p_{T} < 1000 GeV}'%plotlabel

### LOAD PLOTTING UTILITY ###
plot = root.PlotUtility()
plot.Ratio(True) 
plot.SetTDRStyle()
plot.Stack(True)
plot.InitLegend()
plot.SetCut(tcut(cut))
plot.FixRatio(.5)
plot.SetLumi(lumi/1000)
plot.DrawMCErrors(True)
plot.SetNormFactor(False)
plot.AddCMSLabel()
plot.AddLumiLabel()
if plotlabel:
  plot.AddPlotLabel(plotlabel,.18,.77,False,42,.04)
plot.SetMCWeight(weight)

### DEFINE PROCESSES ###
prong1     = root.Process('1-prong',root.kWjets) 
prong2     = root.Process('2-prong',root.kExtra1); prong2.additionalCut = root.TCut("(fj1IsMatched==0||fj1GenSize>1.44)");
prong3     = root.Process('3-prong',root.kTTbar); prong3. additionalCut = root.TCut("(fj1IsMatched==1&&fj1GenSize<1.44)")
qcd        = root.Process('QCD',root.kQCD)
data       = root.Process('Data',root.kData)
data.additionalCut = root.TCut(tOR(sel.triggers['met'],sel.triggers['ele']))
processes = [prong1,prong2,prong3,data]

### ASSIGN FILES TO PROCESSES ###
data.AddFile(basedir+'METSingleElectron.root') 
prong1.AddFile(basedir+'QCD.root')
prong1.AddFile(basedir+'WJets.root')
prong2.AddFile(basedir+'Diboson.root')
prong2.AddFile(basedir+'TTbar.root')
prong2.AddFile(basedir+'SingleTop.root')
prong3.AddFile(basedir+'TTbar.root')
prong3.AddFile(basedir+'SingleTop.root')

#processes = [prong1,data]
for p in processes:
  plot.AddProcess(p)

bins = [40,60,80,100,120,140,160,180,200,220,240,280,350]
nBins = len(bins)-1

msd = root.Distribution("fj1MSD",nBins,"fatjet m_{SD} [GeV]","Events/GeV")
setBins(msd,bins)
#plot.AddDistribution(msd)

plot.AddDistribution(root.Distribution('fj1MSD',50,350,20,'fatjet m_{SD} [GeV]','Events/15 GeV'))
# plot.AddDistribution(root.Distribution('fj1MSDL2L3',50,450,40,'L2L3-corr fatjet m_{SD} [GeV]','Events/10 GeV'))
#plot.AddDistribution(root.Distribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
plot.DrawAll(figsdir+'/'+label)

fhists = root.TFile(figsdir+'/'+label+'hists.root','UPDATE')
h2 = fhists.Get('h_fj1MSD_2-prong')
h1 = fhists.Get('h_fj1MSD_1-prong')
h0 = h2.Clone('h_fj1MSD_bkg')
for iB in xrange(1,h0.GetNbinsX()+1):
  h0.SetBinContent(iB,0)
  h0.SetBinError(iB,0)
h0.Add(h2)
h0.Add(h1)
fhists.WriteTObject(h0)
