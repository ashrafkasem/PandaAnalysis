#!/usr/bin/env python

from os import system,getenv
from sys import argv,exit
import argparse

### SET GLOBAL VARIABLES ###
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--region',metavar='region',type=str,default=None)
args = parser.parse_args()
lumi = 35800.
blind=True
region = args.region
sname = argv[0]

subdirs = {
        'signal' : 'sr',
        'tm' : 'cr_ttbar_mu', 
        'te' : 'cr_ttbar_el', 
        'wmn' : 'cr_w_mu', 
        'wen' : 'cr_w_el', 
        'wmn_fail' : 'cr_w_mu', 
        'wen_fail' : 'cr_w_el', 
        'zmm' : 'cr_dimuon',
        'zmm_fail' : 'cr_dimuon',
        'zee' : 'cr_dielectron',
        'zee_fail' : 'cr_dielectron',
        }

baseDir = getenv('PANDA_FLATDIR')+'/'+subdirs[args.region]+'/' 
dataDir = baseDir#.replace('0_4','0_4_egfix')

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
#import PandaCore.Tools.Functions
import PandaAnalysis.MonoH.Selection_doubleb as sel
#import PandaAnalysis.Monotop.TestSelection as sel
from PandaCore.Drawers.plot_utility import *

### DEFINE REGIONS ###

### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.Stack(True)
plot.Ratio(True)
plot.FixRatio(0.4)
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.cut = sel.cuts[args.region]
plot.SetEvtNum("eventNumber")
plot.SetLumi(lumi/1000)
plot.AddLumiLabel(True)
plot.do_overflow = True
plot.do_underflow = True

weight = sel.weights[region]%lumi
plot.mc_weight = weight

### DEFINE PROCESSES ###
zjets         = Process('Z+jets',root.kZjets)
wjets         = Process('W+jets',root.kWjets)
diboson       = Process('Diboson',root.kDiboson)
ttbar         = Process('t#bar{t}',root.kTTbar)
singletop     = Process('Single t',root.kST)
qcd           = Process("QCD",root.kQCD)
vh            = Process('VH',root.kExtra4)
data          = Process("Data",root.kData)
signal        = Process('m_{V}=1.75 TeV, m_{#chi}=1 GeV',root.kSignal)
processes = [qcd,diboson,singletop,wjets,ttbar,zjets,vh]

### ASSIGN FILES TO PROCESSES ###
if 'signal' in region:
    zjets.add_file(baseDir+'ZtoNuNu.root')
else:
    zjets.add_file(baseDir+'ZJets.root')
wjets.add_file(baseDir+'WJets.root')
diboson.add_file(baseDir+'Diboson.root')
ttbar.add_file(baseDir+'TTbar.root');
singletop.add_file(baseDir+'SingleTop.root')
qcd.add_file(baseDir+'QCD.root')
vh.add_file(baseDir+'VH.root')

if 'e' in region:
    data.additional_cut = sel.triggers['ele']
    data.add_file(dataDir+'SingleElectron.root')
else:
    data.additional_cut = sel.triggers['met']
    data.add_file(dataDir+'MET.root')

if 'signal' not in region:
    processes.append(data)

for p in processes:
    plot.add_process(p)

plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
plot.draw_all(args.outdir+'/'+region+'_')
