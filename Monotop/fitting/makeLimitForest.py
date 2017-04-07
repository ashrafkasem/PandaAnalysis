#!/usr/bin/env python
from re import sub
from sys import argv,exit
from os import path,getenv
from glob import glob
import argparse
parser = argparse.ArgumentParser(description='make forest')
parser.add_argument('--region',metavar='region',type=str,default=None)
toProcess = parser.parse_args().region

argv=[]
import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import *
import PandaCore.Tools.Functions # kinematics
import PandaAnalysis.Monotop.CombinedSelection as sel

Load('PandaAnalysisFlat','LimitTreeBuilder')

baseDir = getenv('PANDA_FLATDIR')
lumi = 36560

factory = root.LimitTreeBuilder()
if toProcess:
  factory.SetOutFile(baseDir+'/limits/limitForest_%s.root'%toProcess)
else:
  factory.SetOutFile(baseDir+'/limits/limitForest_all.root')

def dataCut(basecut,trigger):
  return tAND(trigger,basecut)

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

def shiftBtags(label,tree,varmap,cut,baseweight):
  ps = []
  for shift in ['BUp','BDown','MUp','MDown']:
    for cent in ['sf_btag','sf_sjbtag']:
      shiftedlabel = '_'
      if 'sj' in cent:
        shiftedlabel += 'sj'
      if 'B' in shift:
        shiftedlabel += 'btag'
      else:
        shiftedlabel += 'mistag'
      if 'Up' in shift:
        shiftedlabel += 'Up'
      else:
        shiftedlabel += 'Down'
      weight = sel.weights[baseweight+'_'+cent+shift]%lumi
      shiftedProcess = root.Process(label,tree,varmap,cut,weight)
      shiftedProcess.syst = shiftedlabel
      ps.append(shiftedProcess)
  return ps

# input
tZll,fZll = getTree('ZJets')
tZvv,fZvv = getTree('ZtoNuNu')
tWlv,fWlv = getTree('WJets')
tPho,fPho = getTree('GJets')
tTTbar,fTT = getTree('TTbar')
tVV,fVV = getTree('Diboson')
tQCD,fQCD = getTree('QCD')
tST,fST = getTree('SingleTop')
tMET,fMET = getTree('MET')
tSingleEle,fSEle = getTree('SingleElectron')
tSinglePho,fSPho = getTree('SinglePhoton')
#tSig,fSig = getTree('monotop-nr-v3-1700-100_med-1700_dm-100') # this is just a sample point

tAllSig = {}; fAllSig = {}
if enable('signal'):
  signalFiles = glob(baseDir+'/Vector*root')
  for f in signalFiles:
    fname = f.split('/')[-1].replace('.root','')
    signame = fname
    replacements = {
      'Vector_MonoTop_NLO_Mphi-':'',
      '_gSM-0p25_gDM-1p0_13TeV-madgraph':'',
      '_Mchi-':'_',
      }
    for k,v in replacements.iteritems():
      signame = signame.replace(k,v)
    tAllSig[signame],fAllSig[signame] = getTree(fname)

factory.cd()
regions = {}
processes = {}

vms = {}
for region_type,met_type,phi_type in [('signal','pfmet','pfmetphi'),
                                      ('w','pfUWmag','pfUWphi'),
                                      ('z','pfUZmag','pfUZphi'),
                                      ('a','pfUAmag','pfUAphi')]:
  vms[region_type] = root.VariableMap()
  vms[region_type].AddVar('met',met_type)
  # vms[region_type].AddVar('metphi',phi_type)
  vms[region_type].AddVar('genBosonPt','genBosonPt')
  # vms[region_type].AddVar('genBosonPhi','genBosonPhi')
  for x in ['fj1Tau32','top_ecf_bdt']:
    vms[region_type].AddVar(x,x)


# test region
if enable('test'):
  regions['test'] = root.Region('test')
  cut = sel.cuts['signal']
  weight = sel.weights['signal']%lumi
  processes['test'] = [
    root.Process('Data',tMET,vms['signal'],dataCut(cut,sel.triggers['met']),'1'),
    root.Process('Diboson',tVV,vms['signal'],cut,weight),
  ]
  btag_shifts = []
  for p in processes['test']:
    if p.name=='Data':
      continue
    btag_shifts += shiftCSV(p.name,p.GetInput(),vms['signal'],cut,'signal')
  processes['test'] += btag_shifts
  for p in processes['test']:
    regions['test'].AddProcess(p)
  factory.AddRegion(regions['test'])

# signal region
if enable('signal'):
  regions['signal'] = root.Region('signal')
  cut = sel.cuts['signal']
  weight = sel.weights['signal']%lumi
  vm = vms['signal']
  processes['signal'] = [
    root.Process('Data',tMET,vm,dataCut(cut,sel.triggers['met']),'1'),
    root.Process('Zvv',tZvv,vm,cut,weight),
    root.Process('Zll',tZll,vm,cut,weight),
    root.Process('Wlv',tWlv,vm,cut,weight),
    root.Process('ttbar',tTTbar,vm,cut,weight),
    root.Process('ST',tST,vm,cut,weight),
    root.Process('Diboson',tVV,vm,cut,weight),
    root.Process('QCD',tQCD,vm,cut,weight),
#    root.Process('signal',tSig,vm,cut,weight),
  ]
  for signame,tsig in tAllSig.iteritems():
    processes['signal'].append( root.Process(signame,tsig,vm,cut,weight) )
  btag_shifts = []
  for p in processes['signal']:
    if p.name=='Data':
      continue
    btag_shifts += shiftBtags(p.name,p.GetInput(),vm,cut,'signal')
    #btag_shifts += shiftCSV(p.name,p.GetInput(),vm,cut,'signal')
  processes['signal'] += btag_shifts

  for p in processes['signal']:
    regions['signal'].AddProcess(p)
  factory.AddRegion(regions['signal'])

#singlemuonw
if enable('singlemuonw'):
  regions['singlemuonw'] = root.Region('singlemuonw')
  cut = sel.cuts['singlemuonw']
  weight = sel.weights['singlemuonw']%lumi
  vm = vms['w']
  processes['singlemuonw'] = [
    root.Process('Data',tMET,vm,dataCut(cut,sel.triggers['met']),'1'),
    root.Process('Zvv',tZvv,vm,cut,weight),
    root.Process('Zll',tZll,vm,cut,weight),
    root.Process('Wlv',tWlv,vm,cut,weight),
    root.Process('ttbar',tTTbar,vm,cut,weight),
    root.Process('ST',tST,vm,cut,weight),
    root.Process('Diboson',tVV,vm,cut,weight),
    root.Process('QCD',tQCD,vm,cut,weight),
  ]
  btag_shifts = []
  for p in processes['singlemuonw']:
    if p.name=='Data':
      continue
    btag_shifts += shiftBtags(p.name,p.GetInput(),vm,cut,'singlemuonw')
    #btag_shifts += shiftCSV(p.name,p.GetInput(),vm,cut,'singlemuonw')
  processes['singlemuonw'] += btag_shifts

  for p in processes['singlemuonw']:
    regions['singlemuonw'].AddProcess(p)
  factory.AddRegion(regions['singlemuonw'])

#singleelectronw
if enable('singleelectronw'):
  regions['singleelectronw'] = root.Region('singleelectronw')
  cut = sel.cuts['singleelectronw']
  weight = sel.weights['singleelectronw']%lumi
  vm = vms['w']
  processes['singleelectronw'] = [
    root.Process('Data',tSingleEle,vm,dataCut(cut,sel.triggers['ele']),'1'),
    root.Process('Zvv',tZvv,vm,cut,weight),
    root.Process('Zll',tZll,vm,cut,weight),
    root.Process('Wlv',tWlv,vm,cut,weight),
    root.Process('ttbar',tTTbar,vm,cut,weight),
    root.Process('ST',tST,vm,cut,weight),
    root.Process('Diboson',tVV,vm,cut,weight),
    root.Process('QCD',tQCD,vm,cut,weight),
  ]
  btag_shifts = []
  for p in processes['singleelectronw']:
    if p.name=='Data':
      continue
    btag_shifts += shiftBtags(p.name,p.GetInput(),vm,cut,'singleelectronw')
    #btag_shifts += shiftCSV(p.name,p.GetInput(),vm,cut,'singleelectronw')
  processes['singleelectronw'] += btag_shifts

  for p in processes['singleelectronw']:
    regions['singleelectronw'].AddProcess(p)
  factory.AddRegion(regions['singleelectronw'])

#singlemuontop
if enable('singlemuontop'):
  regions['singlemuontop'] = root.Region('singlemuontop')
  cut = sel.cuts['singlemuontop']
  weight = sel.weights['singlemuontop']%lumi
  vm = vms['w']
  processes['singlemuontop'] = [
    root.Process('Data',tMET,vm,dataCut(cut,sel.triggers['met']),'1'),
    root.Process('Zvv',tZvv,vm,cut,weight),
    root.Process('Zll',tZll,vm,cut,weight),
    root.Process('Wlv',tWlv,vm,cut,weight),
    root.Process('ttbar',tTTbar,vm,cut,weight),
    root.Process('ST',tST,vm,cut,weight),
    root.Process('Diboson',tVV,vm,cut,weight),
    root.Process('QCD',tQCD,vm,cut,weight),
  ]
  btag_shifts = []
  for p in processes['singlemuontop']:
    if p.name=='Data':
      continue
    btag_shifts += shiftBtags(p.name,p.GetInput(),vm,cut,'singlemuontop')
    #btag_shifts += shiftCSV(p.name,p.GetInput(),vm,cut,'singlemuontop')
  processes['singlemuontop'] += btag_shifts

  for p in processes['singlemuontop']:
    regions['singlemuontop'].AddProcess(p)
  factory.AddRegion(regions['singlemuontop'])

#singleelectrontop
if enable('singleelectrontop'):
  regions['singleelectrontop'] = root.Region('singleelectrontop')
  cut = sel.cuts['singleelectrontop']
  weight = sel.weights['singleelectrontop']%lumi
  vm = vms['w']
  processes['singleelectrontop'] = [
    root.Process('Data',tSingleEle,vm,dataCut(cut,sel.triggers['ele']),'1'),
    root.Process('Zvv',tZvv,vm,cut,weight),
    root.Process('Zll',tZll,vm,cut,weight),
    root.Process('Wlv',tWlv,vm,cut,weight),
    root.Process('ttbar',tTTbar,vm,cut,weight),
    root.Process('ST',tST,vm,cut,weight),
    root.Process('Diboson',tVV,vm,cut,weight),
    root.Process('QCD',tQCD,vm,cut,weight),
  ]
  btag_shifts = []
  for p in processes['singleelectrontop']:
    if p.name=='Data':
      continue
    btag_shifts += shiftBtags(p.name,p.GetInput(),vm,cut,'singleelectrontop')
    #btag_shifts += shiftCSV(p.name,p.GetInput(),vm,cut,'singleelectrontop')
  processes['singleelectrontop'] += btag_shifts

  for p in processes['singleelectrontop']:
    regions['singleelectrontop'].AddProcess(p)
  factory.AddRegion(regions['singleelectrontop'])

#dimuon
if enable('dimuon'):
  regions['dimuon'] = root.Region('dimuon')
  cut = sel.cuts['dimuon']
  weight = sel.weights['dimuon']%lumi
  vm = vms['z']
  processes['dimuon'] = [
    root.Process('Data',tMET,vm,dataCut(cut,sel.triggers['met']),'1'),
    root.Process('Zvv',tZvv,vm,cut,weight),
    root.Process('Zll',tZll,vm,cut,weight),
    root.Process('Wlv',tWlv,vm,cut,weight),
    root.Process('ttbar',tTTbar,vm,cut,weight),
    root.Process('ST',tST,vm,cut,weight),
    root.Process('Diboson',tVV,vm,cut,weight),
    root.Process('QCD',tQCD,vm,cut,weight),
  ]

  for p in processes['dimuon']:
    regions['dimuon'].AddProcess(p)
  factory.AddRegion(regions['dimuon'])

#dielectron
if enable('dielectron'):
  regions['dielectron'] = root.Region('dielectron')
  cut = sel.cuts['dielectron']
  weight = sel.weights['dielectron']%lumi
  vm = vms['z']
  processes['dielectron'] = [
    root.Process('Data',tSingleEle,vm,dataCut(cut,sel.triggers['ele']),'1'),
    root.Process('Zvv',tZvv,vm,cut,weight),
    root.Process('Zll',tZll,vm,cut,weight),
    root.Process('Wlv',tWlv,vm,cut,weight),
    root.Process('ttbar',tTTbar,vm,cut,weight),
    root.Process('ST',tST,vm,cut,weight),
    root.Process('Diboson',tVV,vm,cut,weight),
    root.Process('QCD',tQCD,vm,cut,weight),
  ]

  for p in processes['dielectron']:
    regions['dielectron'].AddProcess(p)
  factory.AddRegion(regions['dielectron'])

#photon
if enable('photon'):
  regions['photon'] = root.Region('photon')
  cut = sel.cuts['photon']
  weight = sel.weights['photon']%lumi
  vm = vms['a']
  processes['photon'] = [
    root.Process('Data',tSinglePho,vm,dataCut(cut,sel.triggers['pho']),'1'),
    root.Process('Pho',tPho,vm,cut,weight),
    root.Process('QCD',tSinglePho,vm,dataCut(cut,sel.triggers['pho']),'sf_phoPurity'),
  ]

  for p in processes['photon']:
    regions['photon'].AddProcess(p)
  factory.AddRegion(regions['photon'])

PInfo('makeLimitForest','Starting '+str(toProcess))
factory.Run()
PInfo('makeLimitForest','Finishing '+str(toProcess))

factory.Output()

PInfo('makeLimitForest','Outputted '+str(toProcess))
