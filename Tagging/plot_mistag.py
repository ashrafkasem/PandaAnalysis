#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
import argparse

basedir = getenv('PANDA_FLATDIR')
figsdir = basedir+'/figs'

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--indir',metavar='indir',type=str,default=basedir)
parser.add_argument('--outdir',metavar='outdir',type=str,default=figsdir)
parser.add_argument('--cut',metavar='cut',type=str,default=None)
args = parser.parse_args()

figsdir = args.outdir
basedir = args.indir
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
tcut = root.TCut
Load('Drawers','PlotUtility')

### SET GLOBAL VARIABLES ###
lumi = 12918.
logy=False
cut = 'nFatjet==1 && fj1Pt>250 && fj1MaxCSV<0.46 && nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && UWmag>250 && isojetNBtags==0'
if not args.cut:
  label = 'noCut_'
  plotlabel = '40 GeV < m_{SD}'
  cut = tAND(cut,'fj1MSD>40')
elif args.cut=='mass':
  cut = tAND(cut,'fj1MSD>110 && fj1MSD<210')
  label = 'massCut_'
  plotlabel = '110 < m_{SD} < 210 GeV'
elif args.cut=='massW':
  cut = tAND(cut,'fj1MSD>50 && fj1MSD<100')
  label = 'massWCut_'
  plotlabel = '50 < m_{SD} < 100 GeV'

### LOAD PLOTTING UTILITY ###
plot = root.PlotUtility()
plot.Stack(True)
plot.InitLegend()
plot.SetCut(tcut(cut))
plot.Ratio(1) 
plot.FixRatio(.5)
plot.SetLumi(lumi/1000)
plot.DrawMCErrors(True)
plot.SetTDRStyle()
plot.SetNormFactor(False)
plot.AddCMSLabel()
plot.AddLumiLabel()
if plotlabel:
  plot.AddPlotLabel(plotlabel,.18,.77,False,42,.04)

#weight = '%f*normalizedWeight*sf_pu*sf_lep*sf_ewkV*sf_sjbtag0*sf_btag0*sf_tt'%lumi
weight = '%f*normalizedWeight*sf_pu*sf_lep*sf_ewkV*sf_qcdV*sf_sjbtag0*sf_btag0*sf_tt'%lumi
plot.SetMCWeight(weight)

### DEFINE PROCESSES ###
wjets     = root.Process('W+jets',root.kWjets)
diboson   = root.Process('Diboson',root.kDiboson)
ttbar     = root.Process('t#bar{t} [matched]',root.kTTbar); ttbar.additionalCut = root.TCut('(fj1IsMatched==1&&fj1GenSize<1.44)')
ttbarunmatched     = root.Process('t#bar{t} [unmatched]',root.kExtra1); ttbarunmatched.additionalCut = root.TCut('(fj1IsMatched==0||fj1GenSize>1.44)')
singletop = root.Process('Single t',root.kST)
qcd       = root.Process("QCD",root.kQCD)
data      = root.Process("Data",root.kData); data.additionalCut = root.TCut('(trigger&1)!=0')
processes = [diboson,singletop,wjets,ttbarunmatched,ttbar]

### ASSIGN FILES TO PROCESSES ###
#wjets.AddFile(basedir+'WJets_nlo.root')
wjets.AddFile(basedir+'WJets.root')
diboson.AddFile(basedir+'Diboson.root')
ttbar.AddFile(basedir+'TTbar.root')
ttbarunmatched.AddFile(basedir+'TTbar.root')
singletop.AddFile(basedir+'SingleTop.root')
qcd.AddFile(basedir+'QCD.root')

data.AddFile(basedir+'MET.root') 
processes.append(data)

for p in processes:
  plot.AddProcess(p)

#plot.AddDistribution(root.Distribution('UWmag',250,500,20,'W recoil [GeV]','Events'))

plot.AddDistribution(root.Distribution('top_ecfv6fixed_bdt',-1.,1.,20,'Top ECF BDT v6','Events'))

'''
plot.AddDistribution(root.Distribution('top_ecf_bdt_v2',-0.5,.5,20,'Top ECF BDT v2','Events'))

plot.AddDistribution(root.Distribution('fj1Tau32SD',0,1,20,'Groomed #tau_{32}','Events',999,-999,'tau32SD'))

plot.AddDistribution(root.Distribution('fj1Tau32',0,1,20,'#tau_{32}','Events',999,-999,'tau32'))

plot.AddDistribution(root.Distribution('jet1Pt',15,500,20,'leading jet p_{T} [GeV]','Events'))

plot.AddDistribution(root.Distribution('nJet',-0.5,8.5,9,'N_{jet}','Events'))

plot.AddDistribution(root.Distribution('UWmag',250,500,20,'W recoil [GeV]','Events'))

plot.AddDistribution(root.Distribution('puppimet',0,750,20,'MET [GeV]','Events/37.5 GeV'))

plot.AddDistribution(root.Distribution('fj1Pt',250,1000,20,'fatjet p_{T} [GeV]','Events/37.5 GeV'))

plot.AddDistribution(root.Distribution('fj1MSD',40,450,20,'fatjet m_{SD} [GeV]','Events/12.5 GeV'))

plot.AddDistribution(root.Distribution('fj1ECFN_1_4_10/pow(fj1ECFN_1_3_05,2)',0.6,1.5,20,'_{1}e_{4}^{1.}/(_{1}e_{3}^{0.5})^{2}','Events',999,-999,'input0'))

plot.AddDistribution(root.Distribution('fj1ECFN_2_4_20/pow(fj1ECFN_1_3_20,2)',0.25,3.,20,'_{2}e_{4}^{2.}/(_{1}e_{3}^{02.})^{2}','Events',999,-999,'input1'))

plot.AddDistribution(root.Distribution('fj1ECFN_2_4_10/pow(fj1ECFN_1_3_10,2)',1.,3.,20,'_{2}e_{4}^{1.}/(_{1}e_{3}^{01.})^{2}','Events',999,-999,'input2'))

plot.AddDistribution(root.Distribution('fj1ECFN_3_3_10/pow(fj1ECFN_1_2_20,1.5)',0.,.4,20,'_{3}e_{3}^{1.}/(_{1}e_{2}^{2.})^{3/2}','Events',999,-999,'input3'))

plot.AddDistribution(root.Distribution('fj1ECFN_2_4_10/fj1ECFN_1_2_20',0.,0.015,20,'_{2}e_{4}^{1.}/_{1}e_{2}^{2.}','Events',999,-999,'input4'))

plot.AddDistribution(root.Distribution('fj1ECFN_1_4_20/pow(fj1ECFN_1_3_10,2)',0.25,1.5,20,'_{1}e_{4}^{2.}/(_{1}e_{3}^{1.})^{2}','Events',999,-999,'input5'))

plot.AddDistribution(root.Distribution('fj1ECFN_2_4_05/pow(fj1ECFN_1_3_05,2)',1.25,2.5,20,'_{2}e_{4}^{0.5}/(_{1}e_{3}^{0.5})^{2}','Events',999,-999,'input6'))

plot.AddDistribution(root.Distribution('fj1ECFN_1_3_10/fj1ECFN_2_3_05',0.25,.9,20,'_{1}e_{3}^{1.}/_{2}e_{3}^{0.5}','Events',999,-999,'input7'))

plot.AddDistribution(root.Distribution('fj1ECFN_3_3_10/pow(fj1ECFN_3_3_20,.5)',0.,.35,20,'_{3}e_{3}^{1.}/(_{3}e_{3}^{2.})^{1/2}','Events',999,-999,'input8'))

plot.AddDistribution(root.Distribution('fj1ECFN_3_3_05/pow(fj1ECFN_1_2_05,3.)',1.,3.,20,'_{3}e_{3}^{0.5}/(_{1}e_{2}^{0.5})^{3}','Events',999,-999,'input9'))
'''
'''
plot.AddDistribution(root.Distribution('fj1ECFN_2_4_20/pow(fj1ECFN_1_3_20,2)',0.25,2.5,20,'N_{3}(#beta=2.0)','Events',999,-999,'N3_20'))

plot.AddDistribution(root.Distribution('fj1ECFN_2_4_10/pow(fj1ECFN_1_3_10,2)',0.75,2.5,20,'N_{3}(#beta=1.0)','Events',999,-999,'N3_10'))

plot.AddDistribution(root.Distribution('fj1ECFN_2_4_40/pow(fj1ECFN_1_3_40,2)',0,2,20,'N_{3}(#beta=4.0)','Events',999,-999,'N3_40'))
'''

# plot.AddDistribution(root.Distribution('fj1Tau32SD',0,1,20,'Groomed #tau_{32}','Events',999,-999,'tau32SD'))

# plot.AddDistribution(root.Distribution('fj1Tau32',0,1,20,'#tau_{32}','Events',999,-999,'tau32'))

### DRAW AND CATALOGUE ###
plot.DrawAll(figsdir+'/'+label)
#plot.DrawAll(figsdir+'/nlo_'+label)
