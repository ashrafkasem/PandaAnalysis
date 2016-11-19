#!/usr/bin/env python
from re import sub
from sys import argv,exit
from os import path,getenv
import argparse
parser = argparse.ArgumentParser(description='make forest')
parser.add_argument('--region',metavar='region',type=str,default=None)
toProcess = parser.parse_args().region

argv=[]
import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import *
import PandaCore.Tools.Functions # kinematics
#import PandaAnalysis.VBF.Selection as sel
#import PandaAnalysis.VBF.MonojetSelection as sel
import PandaAnalysis.VBF.LooseSelection as sel

Load('PandaAnalysisFlat','LimitTreeBuilder')

baseDir = getenv('PANDA_ZEYNEPDIR')+'/merged/'
lumi = 12918


factory = root.LimitTreeBuilder()
if toProcess:
  factory.SetOutFile(baseDir+'/limits/limitForest_%s.root'%toProcess)
else:
  factory.SetOutFile(baseDir+'/limits/limitForest_all.root')

def dataCut(basecut,trigger):
  # return tAND('metFilter==1',tAND(trigger,basecut))
  #return tAND(trigger,basecut)
  return tAND(tAND(trigger,basecut),'runNum<=276811')

treelist = []

def getTree(fpath):
  global treelist
  fIn = root.TFile(baseDir+fpath+'.root')
  tIn = fIn.Get('events')
  treelist.append(tIn)
  return tIn,fIn

def enable(regionName):
  if toProcess:
    return (toProcess==regionName)
  else:
    return True

# input
tZll,fZll = getTree('ZJets')
tZvv,fZvv = getTree('ZtoNuNu')
tWlv,fWlv = getTree('WJets')
tWlv_nlo,fWlv_nlo = getTree('WJets_nlo')
tewkZll,fewkZll = getTree('EWKZJets')
tewkZvv,fewkZvv = getTree('EWKZtoNuNu')
tewkWlv,fewkWlv = getTree('EWKWJets')
tPho,fPho = getTree('GJets')
tTTbar,fTT = getTree('TTbar')
tVV,fVV = getTree('Diboson')
tQCD,fQCD = getTree('QCD')
tST,fST = getTree('SingleTop')
tMET,fMET = getTree('MET')
tSingleEle,fSEle = getTree('SingleElectron')
tSinglePho,fSPho = getTree('SinglePhoton')
tVBF,fVBF = getTree('VBF_H125')
tGGF,fGGF = getTree('GGF_H125')

tAllWlv = root.TChain('events')
for f in ['WJets','EWKWJets']:
  tAllWlv.AddFile(baseDir+'/'+f+'.root')
tAllZll = root.TChain('events')
for f in ['ZJets','EWKZJets']:
  tAllZll.AddFile(baseDir+'/'+f+'.root')
tAllZvv = root.TChain('events')
for f in ['ZtoNuNu','EWKZtoNuNu']:
  tAllZvv.AddFile(baseDir+'/'+f+'.root')
treelist += [tAllWlv,tAllZll,tAllZvv]

factory.cd()
regions = {}
processes = {}

vm = root.VariableMap()
vm.AddVar('met','met')
vm.AddVar('metPhi','metPhi')
vm.AddVar('genBosonPt','genBos_pt')
vm.AddVar('genBosonPhi','genBos_phi')
for x in ['jjDEta','mjj','jot1Pt','jot2Pt','jot1Eta','jot2Eta','minJetMetDPhi_withendcap']:
  vm.AddVar(x,x)
vm.AddFormula('jjDPhi','fabs(SignedDeltaPhi(jot1Phi,jot2Phi))')

# test region
if enable('test'):
  regions['test'] = root.Region('test')
  cut = sel.cuts['signal']
  weight = '%f*%s'%(lumi,sel.weights['signal'])
  processes['test'] = [
    root.Process('Data',tMET,vm,dataCut(cut,sel.triggers['met']),'1'),
    root.Process('Diboson',tVV,vm,cut,weight),
  ]

  for p in processes['test']:
    regions['test'].AddProcess(p)
  factory.AddRegion(regions['test'])

# signal region
if enable('signal'):
  regions['signal'] = root.Region('signal')
  cut = sel.cuts['signal']
  weight = '%f*%s'%(lumi,sel.weights['signal'])
  PInfo('makeLimitForest.py',cut)
  processes['signal'] = [
    root.Process('Data',tMET,vm,dataCut(cut,sel.triggers['met']),'1'),
    root.Process('Zvv',tZvv,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('Wlv',tWlv,vm,cut,tTIMES('wkfactor*ewk_w',weight)),
    root.Process('Zll',tZll,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('ewkZvv',tewkZvv,vm,cut,weight),
    root.Process('ewkWlv',tewkWlv,vm,cut,weight),
    root.Process('ewkZll',tewkZll,vm,cut,weight),
    root.Process('ttbar',tTTbar,vm,cut,weight),
    root.Process('ST',tST,vm,cut,weight),
    root.Process('Diboson',tVV,vm,cut,weight),
    root.Process('QCD',tQCD,vm,cut,weight),
    root.Process('VBF_H125',tVBF,vm,cut,weight),
    root.Process('GGF_H125',tGGF,vm,cut,weight),
    root.Process('allWlv',tAllWlv,vm,cut,tTIMES('wkfactor*ewk_w',weight)),
    root.Process('allZvv',tAllZvv,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('allZll',tAllZll,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
  ]

  for p in processes['signal']:
    regions['signal'].AddProcess(p)
  factory.AddRegion(regions['signal'])

# wmn
if enable('wmn'):
  regions['wmn'] = root.Region('wmn')
  cut = sel.cuts['wmn']
  weight = '%f*%s'%(lumi,sel.weights['wmn'])
  processes['wmn'] = [
    root.Process('Data',tMET,vm,dataCut(cut,sel.triggers['met']),'1'),
    root.Process('Zvv',tZvv,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('Wlv',tWlv,vm,cut,tTIMES('wkfactor*ewk_w',weight)),
    root.Process('Zll',tZll,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('ewkZvv',tewkZvv,vm,cut,weight),
    root.Process('ewkWlv',tewkWlv,vm,cut,weight),
    root.Process('ewkZll',tewkZll,vm,cut,weight),
    root.Process('ttbar',tTTbar,vm,cut,weight),
    root.Process('ST',tST,vm,cut,weight),
    root.Process('Diboson',tVV,vm,cut,weight),
    root.Process('QCD',tQCD,vm,cut,weight),
    root.Process('allWlv',tAllWlv,vm,cut,tTIMES('wkfactor*ewk_w',weight)),
    root.Process('allZvv',tAllZvv,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('allZll',tAllZll,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
  ]

  for p in processes['wmn']:
    regions['wmn'].AddProcess(p)
  factory.AddRegion(regions['wmn'])

# wen
if enable('wen'):
  regions['wen'] = root.Region('wen')
  cut = sel.cuts['wen']
  weight = '%f*%s'%(lumi,sel.weights['wen'])
  processes['wen'] = [
    root.Process('Data',tSingleEle,vm,dataCut(cut,sel.triggers['ele']),'1'),
    root.Process('Zvv',tZvv,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('Wlv',tWlv,vm,cut,tTIMES('wkfactor*ewk_w',weight)),
    root.Process('Zll',tZll,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('ewkZvv',tewkZvv,vm,cut,weight),
    root.Process('ewkWlv',tewkWlv,vm,cut,weight),
    root.Process('ewkZll',tewkZll,vm,cut,weight),
    root.Process('ttbar',tTTbar,vm,cut,weight),
    root.Process('ST',tST,vm,cut,weight),
    root.Process('Diboson',tVV,vm,cut,weight),
    root.Process('QCD',tQCD,vm,cut,weight),
    root.Process('allWlv',tAllWlv,vm,cut,tTIMES('wkfactor*ewk_w',weight)),
    root.Process('allZvv',tAllZvv,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('allZll',tAllZll,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
  ]

  for p in processes['wen']:
    regions['wen'].AddProcess(p)
  factory.AddRegion(regions['wen'])

# zmm
if enable('zmm'):
  regions['zmm'] = root.Region('zmm')
  cut = sel.cuts['zmm']
  weight = '%f*%s'%(lumi,sel.weights['zmm'])
  processes['zmm'] = [
    root.Process('Data',tMET,vm,dataCut(cut,sel.triggers['met']),'1'),
    root.Process('Zvv',tZvv,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('Wlv',tWlv,vm,cut,tTIMES('wkfactor*ewk_w',weight)),
    root.Process('Zll',tZll,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('ewkZvv',tewkZvv,vm,cut,weight),
    root.Process('ewkWlv',tewkWlv,vm,cut,weight),
    root.Process('ewkZll',tewkZll,vm,cut,weight),
    root.Process('ttbar',tTTbar,vm,cut,weight),
    root.Process('ST',tST,vm,cut,weight),
    root.Process('Diboson',tVV,vm,cut,weight),
    root.Process('QCD',tQCD,vm,cut,weight),
    root.Process('allWlv',tAllWlv,vm,cut,tTIMES('wkfactor*ewk_w',weight)),
    root.Process('allZvv',tAllZvv,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('allZll',tAllZll,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
  ]

  for p in processes['zmm']:
    regions['zmm'].AddProcess(p)
  factory.AddRegion(regions['zmm'])

# zee
if enable('zee'):
  regions['zee'] = root.Region('zee')
  cut = sel.cuts['zee']
  weight = '%f*%s'%(lumi,sel.weights['zee'])
  processes['zee'] = [
    root.Process('Data',tSingleEle,vm,dataCut(cut,sel.triggers['ele']),'1'),
    root.Process('Zvv',tZvv,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('Wlv',tWlv,vm,cut,tTIMES('wkfactor*ewk_w',weight)),
    root.Process('Zll',tZll,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('ewkZvv',tewkZvv,vm,cut,weight),
    root.Process('ewkWlv',tewkWlv,vm,cut,weight),
    root.Process('ewkZll',tewkZll,vm,cut,weight),
    root.Process('ttbar',tTTbar,vm,cut,weight),
    root.Process('ST',tST,vm,cut,weight),
    root.Process('Diboson',tVV,vm,cut,weight),
    root.Process('QCD',tQCD,vm,cut,weight),
    root.Process('allWlv',tAllWlv,vm,cut,tTIMES('wkfactor*ewk_w',weight)),
    root.Process('allZvv',tAllZvv,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
    root.Process('allZll',tAllZll,vm,cut,tTIMES('zkfactor*ewk_z',weight)),
  ]

  for p in processes['zee']:
    regions['zee'].AddProcess(p)
  factory.AddRegion(regions['zee'])


# photon
if enable('pho'):
  regions['pho'] = root.Region('pho')
  cut = sel.cuts['pho']
  weight = '%f*%s'%(lumi,sel.weights['pho'])
  processes['pho'] = [
    root.Process('Data',tSinglePho,vm,dataCut(cut,sel.triggers['pho']),'1'),
    root.Process('Pho',tPho,vm,cut,tTIMES('akfactor*ewk_a',weight)),
    # root.Process('QCD',tSinglePho,vmA,dataCut(cut,phoTrigger),'photonPurityWeight'),
    root.Process('QCD',tQCD,vm,cut,weight),
  ]
  for p in processes['pho']:
    regions['pho'].AddProcess(p)
  factory.AddRegion(regions['pho'])


PInfo('makeLimitForest','Starting '+str(toProcess))
factory.Run()
PInfo('makeLimitForest','Finishing '+str(toProcess))

for t in treelist:
  t.SetDirectory(0)

factory.Output()

