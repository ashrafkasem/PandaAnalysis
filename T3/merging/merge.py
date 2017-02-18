#!/usr/bin/env python

from array import array
from glob import glob
from re import sub
from sys import argv
from os import environ,system,path

sname = argv[0]
arguments = [x for x in argv[1:]] # deep copy
argv=[]

from ROOT import gSystem, gROOT
import ROOT as root
from PandaCore.Tools.process import *
from PandaCore.Tools.Misc import *

gROOT.LoadMacro("${CMSSW_BASE}/src/PandaCore/Tools/interface/Normalizer.h")
gSystem.Load('libPandaCoreTools.so')

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

inbase = environ['SUBMIT_OUTDIR']
outbase = environ['PANDA_FLATDIR']

#suffix = ' > /dev/null '
suffix = ''

def hadd(inpath,outpath):
  if type(inpath)==type('str'):
    infiles = glob(inpath)
    PInfo(sname,'hadding %s into %s'%(inpath,outpath))
    cmd = 'hadd -k -ff -n 100 -f %s %s %s'%(outpath,inpath,suffix)
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
  system(cmd+suffix)

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
  PInfo(sname,'normalizing %s (%s) ...'%(fpath,opt))
  n = root.Normalizer();
  n.NormalizeTree(fpath,xsec)

def merge(shortnames,mergedname):
  for shortname in shortnames:
    if 'monotop' in shortname:
      pd = shortname
      xsec = 1
    elif shortname in pds:
      pd = pds[shortname][0]
      xsec = pds[shortname][1]
    else:
      for shortname_ in [shortname.split('_')[0],shortname.split('_')[-1]]:
        if shortname_ in pds:
          pd = pds[shortname_][0]
          xsec = pds[shortname_][1]
          break
    inpath = inbase+shortname+'_*.root'
    hadd(inpath,'/tmp/%s/split/%s.root'%(user,shortname))
    if xsec>0:
      normalizeFast('/tmp/%s/split/%s.root'%(user,shortname),xsec)
  hadd(['/tmp/%s/split/%s.root'%(user,x) for x in shortnames],'/tmp/%s/merged/%s.root'%(user,mergedname))

d = {
  'test'                : ['Diboson_ww'],
  'Diboson'             : ['Diboson_ww','Diboson_wz','Diboson_zz'],
  'ZJets'               : ['ZJets_ht%sto%s'%(str(x[0]),str(x[1])) for x in [(100,200),(200,400),(400,600),(600,800),(800,1200),(1200,2500),(2500,'inf')]],
  'ZtoNuNu'             : ['ZtoNuNu_ht100to200','ZtoNuNu_ht200to400','ZtoNuNu_ht400to600','ZtoNuNu_ht600to800','ZtoNuNu_ht800to1200','ZtoNuNu_ht1200to2500','ZtoNuNu_ht2500toinf'],
  'GJets'               : ['GJets_ht100to200','GJets_ht200to400','GJets_ht400to600','GJets_ht600toinf'],
  'WJets'               : ['WJets_ht100to200','WJets_ht200to400','WJets_ht400to600','WJets_ht600to800','WJets_ht800to1200','WJets_ht1200to2500','WJets_ht2500toinf'],
  'TTbar'               : ['TTbar_Powheg'],
  'TTbar_isrup'         : ['TTbar_PowhegISRUp'],
  'TTbar_isrdown'       : ['TTbar_PowhegISRDown'],
  'TTbar_tuneup'        : ['TTbar_PowhegTuneUp'],
  'TTbar_tunedown'      : ['TTbar_PowhegTuneDown'],
  'TTbar_FXFX'		: ['TTbar_FXFX'],
  'TTbar_Herwig'        : ['TTbar_Herwig'],
  'SingleTop'           : ['SingleTop_tT','SingleTop_tTbar','SingleTop_tbarW','SingleTop_tW'],
  'QCD'                 : ['QCD_ht100to200','QCD_ht200to300','QCD_ht300to500','QCD_ht500to700','QCD_ht700to1000','QCD_ht1000to1500','QCD_ht1500to2000','QCD_ht2000toinf'],
  'MET'                 : ['MET'],
  'SingleElectron'      : ['SingleElectron'],
  'DoubleEG'            : ['DoubleEG'],
  'SinglePhoton'        : ['SinglePhoton'],
  'ZJets_nlo'           : ['ZJets_nlo'],
  'WJets_nlo'           : ['WJets_pt%sto%s'%(str(x[0]),str(x[1])) for x in [(100,250),(250,400),(400,600),(600,'inf')] ],
  'ZHbb'                : ['ZHbb_mH125'],
  'ggZHbb'              : ['ggZHbb_mH125'],
  'WpH'                 : ['WpLNuHbb'],
  'WmH'                 : ['WmLNuHbb'],
}

args = {}

for pd in arguments:
  if pd in d:
    args[pd] = d[pd]
  else:
    args[pd] = [pd]

for pd in args:
  merge(args[pd],pd)
  system('cp -r /tmp/%s/merged/%s.root %s'%(user,pd,outbase))
  PInfo(sname,'finished with '+pd)

