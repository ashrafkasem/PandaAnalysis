#!/usr/bin/env python

from array import array
from glob import glob
from re import sub
from sys import argv
from os import environ,system,path

sname = argv[0]
arguments = [x for x in argv[1:]]
argv=[]

import ROOT as root
from PandaCore.Tools.process import *
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import *

Load('Tools','Normalizer')

pds = {}
for k,v in processes.iteritems():
  if v[1]=='MC':
    pds[v[0]] = (k,v[2])  
  else:
    pds[v[0]] = (k,-1)

VERBOSE=False

user = environ['USER']
system('mkdir -p /tmp/%s/split'%user) # tmp dir
system('mkdir -p /tmp/%s/merged'%user) # tmp dir

inbase = environ['PANDA_ZEYNEPDIR']
#inbase = environ['PANDA_ZEYNEPDIR_PROMPT']
outbase = inbase + '/merged/'

def hadd(inpath,outpath):
  if type(inpath)==type('str'):
    infiles = glob(inpath)
    PInfo(sname,'hadding %s into %s'%(inpath,outpath))
    cmd = 'hadd -k -ff -n 100 -f %s %s > /dev/null'%(outpath,inpath)
    system(cmd)
    return
  else:
    infiles = inpath
  if len(infiles)==0:
    PWarning(sname,'nothing hadded into',outpath)
    return
  elif len(infiles)==1:
    cmd = 'cp %s %s'%(infiles[0],outpath)
  else:
    cmd = 'hadd -k -ff -n 100 -f %s '%outpath
    for f in infiles:
      if path.isfile(f):
        cmd += '%s '%f
  if VERBOSE: PInfo(sname,cmd)
  system(cmd+' >/dev/null 2>/dev/null')

def normalizeFast(fpath,opt):
  xsec=-1
  if type(opt)==type(1.) or type(opt)==type(1):
    xsec = opt
  else:
    try:
      xsec = processes[proc][2]
    except KeyError:
      for k,v in processes.iteritems():
        if proc in k:
          xsec = v[2]
  if xsec<0:
    PError(sname,'could not find xsec, skipping %s!'%opt)
    return
  PInfo(sname,'normalizing %s (%s) ...'%(fpath,xsec))
  n = root.Normalizer()
  n.histName = 'htotal'
  print opt,xsec
  n.NormalizeTree(fpath,xsec)

def merge(shortnames,mergedname):
  for shortname in shortnames:
    myname=None
    for k,v in processes.iteritems():
      if 'Zprime' in k:
        continue  
      if shortname in k:
        myname = v[0]
        isData = (v[1]=='Data')
        xsec = v[2]
        break
    if myname==None:
      PError(sname,'Could not translate %s'%shortname)
      continue
    if not isData:
      normalizeFast(inbase+'monojet_'+shortname+'.root',xsec)
  hadd([inbase+'monojet_'+x+'.root' for x in shortnames],'/tmp/%s/merged/%s.root'%(user,mergedname))


d = {
  'Diboson'             : ['WW','WZ','ZZ'],
  'VBF_H125'            : ['VBF_HToInvisible_M125_13TeV_powheg_pythia8'],
  'GGF_H125'            : ['Glu_HToInvisible_M125_13TeV_powheg_pythia8'],
  'ZJets'               : ['DYJetsToLL_M-50_HT-%sto%s'%(str(x[0]),str(x[1])) for x in [(100,200),(200,400),(400,600),(600,800),(800,1200),(1200,2500),(2500,'Inf')]],
  'ZtoNuNu'             : ['ZJetsToNuNu_HT-%sTo%s_13TeV'%(str(x[0]),str(x[1])) for x in [(100,200),(200,400),(400,600),(600,800),(800,1200),(1200,2500),(2500,'Inf')]],
  'WJets'               : ['WJetsToLNu_HT-%sTo%s'%(str(x[0]),str(x[1])) for x in [(100,200),(200,400),(400,600),(600,800),(800,1200),(1200,2500),(2500,'Inf')]],
  'WJets_nlo'           : ['WJetsToLNu_Pt-%sTo%s'%(str(x[0]),str(x[1])) for x in [(100,250),(250,400),(400,600),(600,'Inf')]],
  'TTbar'               : ['TT'],
  'MET'                 : ['MET'],
  'SingleElectron'      : ['SingleElectron'],
  'SinglePhoton'        : ['SinglePhoton'],
  'EWKZtoNuNu'          : ['EWKZ2Jets_ZToNuNu_13TeV'],
  'EWKZJets'            : ['EWKZ2Jets_ZToLL_M-50_13TeV'],
  'EWKWJets'            : ['EWKW%s2Jets_WToLNu_M-50_13TeV'%(x) for x in ['Plus','Minus']],
  'SingleTop'           : ['ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8','ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8','ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8','ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8'],
  'QCD'                 : ['QCD_HT%sto%s'%(str(x[0]),str(x[1])) for x in [(200,300),(300,500),(500,700),(700,1000),(1000,1500),(1500,2000),(2000,'Inf')]],
  'GJets'               : ['GJets_HT-%sTo%s'%(str(x[0]),str(x[1])) for x in [(100,200),(200,400),(400,600),(600,'Inf')]],
}

args = {}

for pd in arguments:
  args[pd] = d[pd]

for pd in args:
  merge(args[pd],pd)
  system('cp /tmp/%s/merged/%s.root %s'%(user,pd,outbase))
  PInfo(sname,'finished with '+pd)

