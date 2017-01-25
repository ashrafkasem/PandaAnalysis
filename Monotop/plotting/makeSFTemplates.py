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
parser.add_argument('--jesr',metavar='jesr',type=str,default='central')
parser.add_argument('--cut',metavar='cut',type=float,default=0.5)
args = parser.parse_args()

figsdir = args.outdir + '/' + args.jesr
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

cut = sel.cuts['singlemuontop']
#cut = tOR(sel.cuts['singlemuontop'],sel.cuts['singleelectrontop'])
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
	if 'Tau' in args.disc:
		cut = tAND(cut,'%s<%f'%(args.disc,args.cut))
		plotlabel = '%s<%.2f'%(disclabels[args.disc],args.cut)
	else:
		cut = tAND(cut,'%s>%f'%(args.disc,args.cut))
		plotlabel = '%s>%.2f'%(disclabels[args.disc],args.cut)
else:
	label += 'fail_'
	if 'Tau' in args.disc:
		cut = tAND(cut,'%s>%f'%(args.disc,args.cut))
		plotlabel = '%s>%.2f'%(disclabels[args.disc],args.cut)
	else:
		cut = tAND(cut,'%s<%f'%(args.disc,args.cut))
		plotlabel = '%s<%.2f'%(disclabels[args.disc],args.cut)

cut = tAND(cut,'fj1Pt>250 && fj1Pt<1000')
pt_name = 'fj1Pt'
msd_name = 'fj1MSD'
if args.jesr == 'smeared':
	pt_name = 'fj1PtSmeared'
	msd_name = 'fj1MSDSmeared'
	plotlabel = '#splitline{%s}{Smeared p_{T} and m_{SD}}'%plotlabel
elif args.jesr == 'scaleUp':
	pt_name = 'fj1PtScaleUp'
	msd_name = 'fj1MSDScaleUp'
	plotlabel = '#splitline{%s}{JES Up}'%plotlabel
elif args.jesr == 'scaleDown':
	pt_name = 'fj1PtScaleDown'
	msd_name = 'fj1MSDScaleDown'
	plotlabel = '#splitline{%s}{JES Down}'%plotlabel
cut = cut.replace('fj1Pt',pt_name)
cut = cut.replace('fj1MSD',msd_name)

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
data.additionalCut = root.TCut(sel.triggers['met'])
processes = [prong1,prong2,prong3,data]

### ASSIGN FILES TO PROCESSES ###
#data.AddFile(basedir+'METSingleElectron.root') 
data.AddFile(basedir+'MET.root') 
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

msd = root.Distribution(msd_name,nBins,"fatjet m_{SD} [GeV]","Events/GeV")
setBins(msd,bins)
#plot.AddDistribution(msd)

plot.AddDistribution(root.Distribution(msd_name,50,350,20,'fatjet m_{SD} [GeV]','Events/15 GeV'))
plot.AddDistribution(root.Distribution(pt_name,250,1000,20,'fatjet p_{T} [GeV]','Events/37.5 GeV'))
#plot.AddDistribution(root.Distribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
plot.DrawAll(figsdir+'/'+label)

fhists = root.TFile(figsdir+'/'+label+'hists.root','UPDATE')
h2 = fhists.Get('h_%s_2-prong'%msd_name)
h1 = fhists.Get('h_%s_1-prong'%msd_name)
h0 = h2.Clone('h_%s_bkg'%msd_name)
for iB in xrange(1,h0.GetNbinsX()+1):
	h0.SetBinContent(iB,0)
	h0.SetBinError(iB,0)
h0.Add(h2)
h0.Add(h1)
fhists.WriteTObject(h0)
