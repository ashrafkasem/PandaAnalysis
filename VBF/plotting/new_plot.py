#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_ZEYNEPDIR')+'/merged/' 
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--basedir',metavar='basedir',type=str,default=None)
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
parser.add_argument('--region',metavar='region',type=str,default=None)
parser.add_argument('--tt',metavar='tt',type=str,default='')
args = parser.parse_args()
lumi = 36560.
blind=True
linear=False
region = args.region
sname = argv[0]
if args.basedir:
	baseDir = args.basedir

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.VBF.LooseSelection as sel
from PandaCore.Drawers.plot_utility import *

### DEFINE REGIONS ###

cut = tAND(sel.cuts[args.region],args.cut)

PInfo(sname,'using cut: '+cut)

### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.Stack(True)
plot.Logy(not(linear))
if ('signal' in region) and blind:
  plot.SetSignalScale(10)
  plot.Ratio(False)
else:
  plot.SetLumi(lumi/1000)
  plot.Ratio(True)
  plot.FixRatio(0.4)
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.cut = cut
plot.AddLumiLabel(True)
plot.do_overflow = True
plot.do_underflow = True

weight = sel.weights[region]%lumi
plot.mc_weight = weight

### DEFINE PROCESSES ###
zjets     = Process('QCD Z+jets',root.kZjets)
wjets     = Process('QCD W+jets',root.kWjets)
zewk     = Process('EWK Z+jets',root.kExtra2)
wewk     = Process('EWK W+jets',root.kExtra3)
diboson   = Process('Diboson',root.kDiboson)
top     = Process('Top',root.kTTbar)
qcd       = Process("QCD",root.kQCD)
data      = Process("Data",root.kData)
processes = [qcd,diboson,top,wewk,zewk,wjets,zjets]

### deal with V+jets ###
zjets.additional_weight = 'zkfactor*ewk_z'
wjets.additional_weight = 'wkfactor*ewk_w'
# hack the old EWK
zewk.use_common_weight = False
zewk.additional_weight = weight.replace('_v2','').replace('_id','').replace('lepton_SF1_iso','1').replace('lepton_SF2_iso','1').replace('EleTrigger','')
wewk.use_common_weight = False
wewk.additional_weight = weight.replace('_v2','').replace('_id','').replace('lepton_SF1_iso','1').replace('lepton_SF2_iso','1').replace('EleTrigger','')

### ASSIGN FILES TO PROCESSES ###
if 'signal' in region:
  zjets.add_file(baseDir+'ZtoNuNu.root')
  zewk.add_file(baseDir+'EWKZtoNuNu.root')
else:
  zjets.add_file(baseDir+'ZJets.root')
  zewk.add_file(baseDir+'EWKZJets.root')
wjets.add_file(baseDir+'WJets.root')
wewk.add_file(baseDir+'EWKWJets.root')
diboson.add_file(baseDir+'Diboson.root')
top.add_file(baseDir+'TTbar.root')
top.add_file(baseDir+'SingleTop.root')
qcd.add_file(baseDir+'QCD.root')

if region in ['signal','wmn','zmm']:
  data.additional_cut = sel.triggers['met']
  data.add_file(baseDir+'MET.root')
  lep='#mu'
elif region in ['zee','wen']:
  if 'z' in region:
    data.additional_cut = tOR(sel.triggers['ele'],sel.triggers['pho'])
  else:
    data.additional_cut = sel.triggers['ele']
  data.add_file(baseDir+'SingleElectron.root')
  lep='e'
if 'signal' not in region:
	processes.append(data)

for p in processes:
  plot.add_process(p)

recoilBins = [200., 230., 260.0, 290.0, 320.0, 350.0, 390.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]

### CHOOSE DISTRIBUTIONS, LABELS ###
if 'signal' in region:
  recoil=VDistribution("met",recoilBins,"MET [GeV]","Events/GeV")
elif region in ['wmn','wen']: 
  plot.add_distribution(FDistribution('lep1Pt',0,500,20,'lep 1 p_{T} [GeV]','Events/25 GeV'))
  plot.add_distribution(FDistribution('lep1Eta',-2.5,2.5,20,'lep 1 #eta','Events'))
  recoil=VDistribution("met",recoilBins,"U(%s) [GeV]"%(lep),"Events/GeV")
elif region in ['zmm','zee']: 
  plot.add_distribution(FDistribution('lep1Pt',0,500,20,'lep 1 p_{T} [GeV]','Events/25 GeV'))
  plot.add_distribution(FDistribution('lep1Eta',-2.5,2.5,20,'lep 1 #eta','Events'))
  plot.add_distribution(FDistribution('lep2Pt',0,500,20,'lep 2 p_{T} [GeV]','Events/25 GeV'))
  plot.add_distribution(FDistribution('lep2Eta',-2.5,2.5,20,'lep 2 #eta','Events'))
  recoil=VDistribution("met",recoilBins,"U(%s%s) [GeV]"%(lep,lep),"Events/GeV")
plot.add_distribution(recoil)
if 'signal' not in region:
	plot.add_distribution(FDistribution('trueMet',0,1000,10,'MET [GeV]','Events/100 GeV'))

plot.add_distribution(FDistribution('mjj',0,4000,20,'m_{jj} [GeV]','Events/200 GeV'))
plot.add_distribution(FDistribution('jjDEta',0,10,20,'#Delta#eta(j_{1},j_{2})','Events'))
plot.add_distribution(FDistribution("fabs(SignedDeltaPhi(jot1Phi,jot2Phi))",0,3.142,20,"#Delta #phi leading jets","Events",filename='jjDPhi'))
plot.add_distribution(FDistribution("jot1Eta",-5,5,20,"Jet 1 #eta","Events"))
plot.add_distribution(FDistribution("jot2Eta",-5,5,20,"Jet 2 #eta","Events"))
plot.add_distribution(FDistribution("jot1Pt",80,500,20,"Jet 1 p_{T} [GeV]","Events"))
plot.add_distribution(FDistribution("jot2Pt",40,500,20,"Jet 2 p_{T} [GeV]","Events"))
plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
plot.draw_all(args.outdir+'/'+region+'_')
