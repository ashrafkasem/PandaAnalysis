#!/usr/bin/env python

from os import system,getenv
from sys import argv,exit
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
parser.add_argument('--region',metavar='region',type=str,default=None)
parser.add_argument('--era',metavar='era',default=None)
parser.add_argument('--tt',metavar='tt',type=str,default=None)
args = parser.parse_args()
lumi = 36560.
blind=True
linear=False
region = args.region
sname = argv[0]
print baseDir

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.Monotop.MonojetSelection as sel
Load('Drawers','PlotUtility')

### DEFINE REGIONS ###

cut = tAND(sel.cuts[args.region],args.cut)
datacut = '1==1'
if args.era:
	all_runs = {
			'B' : (272007,275376),
			'C' : (275657,276283),
			'D' : (276315,276811),
			'E' : (276831,277420),
			'F' : (277772,278808),
			'G' : (278820,280385),
			'H' : (280919,284044),
			}

	run_boundaries = []
	for e in args.era:
		lo,hi = all_runs[e]
		run_boundaries += [lo,hi]
	runs = (min(run_boundaries),max(run_boundaries))
	datacut = 'runNumber>%i && runNumber<%i'%(runs[0],runs[1])
PInfo(sname,'using cut %s'%(cut))

### LOAD PLOTTING UTILITY ###
plot = root.PlotUtility()
plot.Stack(True)
plot.Logy(not(linear))
# plot.SetSignalScale(10)
plot.Ratio(True)
plot.FixRatio(0.4)
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.SetCut(root.TCut(cut))
plot.SetEvtNum("eventNumber")
if ('signal' in region) and blind:
	plot.SetEvtMod(5)
	plot.SetLumi(lumi/5000)
	plot.AddPlotLabel('Every 5th event',.18,.7,False,42,.04)
else:
	plot.SetLumi(lumi/1000)
plot.AddLumiLabel(True)
plot.SetDoOverflow()
plot.SetDoUnderflow()
if args.era:
	plot.SetNormFactor(True)

weight = sel.weights[region]%lumi
plot.SetMCWeight(root.TCut(weight))
PInfo(sname,'using weight: '+weight)


### DEFINE PROCESSES ###
zjets		 = root.Process('Z+jets',root.kZjets)
wjets		 = root.Process('W+jets',root.kWjets)
diboson	 = root.Process('Diboson',root.kDiboson)
ttbar		 = root.Process('t#bar{t}',root.kTTbar); 
if args.tt:
	ttbar.additionalWeight = root.TCut(args.tt)
singletop = root.Process('Single t',root.kST)
qcd			 = root.Process("QCD",root.kQCD)
gjets		 = root.Process('#gamma+jets',root.kGjets)
data			= root.Process("Data",root.kData)
signal		= root.Process('m_{V}=1.7 TeV, m_{#chi}=100 GeV',root.kSignal)
processes = [qcd,diboson,singletop,wjets,ttbar,zjets]

### ASSIGN FILES TO PROCESSES ###
if 'signal' in region:
	zjets.AddFile(baseDir+'ZtoNuNu.root')
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
	qcd.additionalCut = root.TCut(tAND(datacut,sel.triggers['pho']))
	qcd.useCommonWeight = False
	qcd.additionalWeight = root.TCut('sf_phoPurity')
else:
	qcd.AddFile(baseDir+'QCD.root')

if any([x in region for x in ['singlemuonw','singleelectronw']]):
	processes = [qcd,diboson,singletop,zjets,ttbar,wjets,]
if any([x in region for x in ['singlemuontop','singleelectrontop']]):
	processes = [qcd,diboson,singletop,zjets,wjets,ttbar]
if any([x in region for x in ['signal','muon']]):
	data.additionalCut = root.TCut(tAND(datacut,sel.triggers['met']))
	PInfo(sname,'Using MET data')
	data.AddFile(baseDir+'MET.root')
	lep='#mu'
elif 'electron' in region:
	if 'di' in region:
		data.additionalCut = root.TCut(tAND(datacut,tOR(sel.triggers['ele'],sel.triggers['pho'])))
	else:
#		pass
		data.additionalCut = root.TCut(tAND(datacut,sel.triggers['ele']))
	data.AddFile(baseDir+'/SingleElectron.root')
	lep='e'
elif region=='photon':
	data.additionalCut = root.TCut(tAND(datacut,sel.triggers['pho']))
	data.AddFile(baseDir+'SinglePhoton.root')
processes.append(data)

for p in processes:
	plot.AddProcess(p)

#recoilBins = [250,280,310,350,400,450,600,1000]
recoilBins = [200., 230., 260.0, 290.0, 320.0, 350.0, 390.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]
nRecoilBins = len(recoilBins)-1

### CHOOSE DISTRIBUTIONS, LABELS ###
recoil=None

pfmet=root.Distribution("pfmet",nRecoilBins,"PF E_{T}^{miss} [GeV]","Events/GeV")
if 'signal' in region:
	recoil=pfmet
	recoil=root.Distribution("pfmet",nRecoilBins,"PF E_{T}^{miss} [GeV]","Events/GeV")
	dphi=root.Distribution("dphipfmet",0,3.14,20,"min#Delta#phi(jet,PF E_{T}^{miss})","Events")
	plot.AddDistribution(root.Distribution('fabs(calomet-pfmet)/pfmet',0,1,20,'|E_{T,calo}^{miss}-E_{T}^{miss}|/E_{T}^{miss}','Events',999,-999,'pfcalobalance'))
elif any([x in region for x in ['singlemuonw','singleelectronw','singlemuontop','singleelectrontop','singlemuon','singleelectron']]):
	recoil=root.Distribution('pfUWmag',nRecoilBins,'U(%s) [GeV]'%lep,"Events/GeV")
	dphi=root.Distribution("dphipfUW",0,3.14,20,"min#Delta#phi(jet,U(%s))"%lep,"Events")
	plot.AddDistribution(root.Distribution("dphipfmet",0,3.14,20,"min#Delta#phi(jet,PF E_{T}^{miss})","Events"))
	plot.AddDistribution(root.Distribution('looseLep1Pt',0,500,20,'Leading lep p_{T} [GeV]','Events/25 GeV'))
	plot.AddDistribution(root.Distribution('looseLep1Eta',-2.5,2.5,20,'Leading lep #eta','Events'))
	plot.AddDistribution(root.Distribution('fabs(calomet-pfmet)/pfUWmag',0,1,20,'|E_{T,calo}^{miss}-E_{T}^{miss}|/U(%s)'%lep,'Events',999,-999,'pfcalobalance'))
elif any([x in region for x in ['dielectron','dimuon']]):
	recoil=root.Distribution('pfUZmag',nRecoilBins,'U(%s%s) [GeV]'%(lep,lep),"Events/GeV")
	dphi=root.Distribution("dphiUZ",0,3.14,20,"min#Delta#phi(jet,U(%s%s))"%(lep,lep),"Events")
	plot.AddDistribution(root.Distribution("dphipfmet",0,3.14,20,"min#Delta#phi(jet,PF E_{T}^{miss})","Events"))
	plot.AddDistribution(root.Distribution('looseLep1Pt',0,500,20,'Leading lep p_{T} [GeV]','Events/25 GeV'))
	plot.AddDistribution(root.Distribution('looseLep1Eta',-2.5,2.5,20,'Leading lep #eta','Events'))
	plot.AddDistribution(root.Distribution('looseLep2Pt',0,500,20,'Sub-leading lep p_{T} [GeV]','Events/25 GeV'))
	plot.AddDistribution(root.Distribution('looseLep2Eta',-2.5,2.5,20,'Sub-leading lep #eta','Events'))
	plot.AddDistribution(root.Distribution('fabs(calomet-pfmet)/pfUZmag',0,1,20,'|E_{T,calo}^{miss}-E_{T}^{miss}|/U(%s%s)'%(lep,lep),'Events',999,-999,'pfcalobalance'))
elif region=='photon':
	recoil=root.Distribution('pfUAmag',nRecoilBins,'U(#gamma) [GeV]',"Events/GeV")
	dphi=root.Distribution("dphiUA",0,3.14,20,"min#Delta#phi(jet,U(#gamma))","Events")
	plot.AddDistribution(root.Distribution("dphipfmet",0,3.14,20,"min#Delta#phi(jet,pf E_{T}^{miss})","Events"))
	plot.AddDistribution(root.Distribution('loosePho1Pt',175,1175,20,'Leading pho p_{T} [GeV]','Events/50 GeV'))
	plot.AddDistribution(root.Distribution('loosePho1Eta',-2.5,2.5,20,'Leading pho #eta','Events'))
	plot.AddDistribution(root.Distribution('loosePho1Phi',-3.142,3.142,20,'Leading pho #phi','Events'))
	plot.AddDistribution(root.Distribution('fabs(calomet-pfmet)/pfUAmag',0,1,20,'|E_{T,calo}^{miss}-E_{T}^{miss}|/U(#gamma)','Events',999,-999,'pfcalobalance'))
#if region!='photon':
	#plot.AddDistribution(root.Distribution('fixed_mt',0,1000,20,'M_{T} [GeV]','Events/50 GeV'))
if recoil:
	setBins(recoil,recoilBins)
	plot.AddDistribution(recoil)
if recoil!=pfmet:
	setBins(pfmet,recoilBins)
	plot.AddDistribution(pfmet)
plot.AddDistribution(dphi)


plot.AddDistribution(root.Distribution('pfmetphi',-3.142,3.142,20,'PF E_{T}^{miss} #phi','Events'))

plot.AddDistribution(root.Distribution('top_ecf_bdt',-1,1,20,'ECF+#tau_{32}^{SD}+f_{rec} BDT','Events'))

plot.AddDistribution(root.Distribution('fj1Tau32SD',0,1,20,'Groomed #tau_{32}','Events',999,-999,'tau32SD'))

plot.AddDistribution(root.Distribution('fj1Tau32',0,1,20,'#tau_{32}','Events',999,-999,'tau32'))

plot.AddDistribution(root.Distribution('jet1Pt',15,1000,20,'leading jet p_{T} [GeV]','Events'))

plot.AddDistribution(root.Distribution('nJet',-0.5,8.5,9,'N_{jet}','Events'))

plot.AddDistribution(root.Distribution('fj1Pt',250,1000,20,'fatjet p_{T} [GeV]','Events/37.5 GeV'))

plot.AddDistribution(root.Distribution('fj1Eta',-2.5,2.5,20,'fatjet #eta','Events'))

plot.AddDistribution(root.Distribution('fj1MSD',0,250,25,'fatjet m_{SD} [GeV]','Events/10 GeV'))

plot.AddDistribution(root.Distribution('fj1MaxCSV',0,1,20,'fatjet max subjet CSV','Events',999,-999,'fj1MaxCSV'))

plot.AddDistribution(root.Distribution('jet1CSV',0,1,20,'jet 1 CSV','Events',999,-999,'jet1CSV'))


plot.AddDistribution(root.Distribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
if args.era:
	plot.DrawAll(args.outdir+'/'+region+'_'+args.era+'_')
else:
	plot.DrawAll(args.outdir+'/'+region+'_')
