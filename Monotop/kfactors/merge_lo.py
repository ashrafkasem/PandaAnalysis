#!/usr/bin/env python

from array import array
from glob import glob
from re import sub
from sys import argv,exit
from os import environ,system,path

sname = argv[0]
arguments = [x for x in argv[1:]] # deep copy
argv=[]

import ROOT as root
from PandaCore.Tools.process import *
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import Load

Load('Normalizer')

translate = {
        'GJets_HT-100to200':'GJets_ht100to200',
        'GJets_HT-200to400':'GJets_ht200to400',
        'GJets_HT-400to600':'GJets_ht400to600',
        'GJets_HT-600toInf':'GJets_ht600toinf',

        'DYJets_HT-100to200':'ZJets_ht100to200',
        'DYJets_HT-200to400':'ZJets_ht200to400',
        'DYJets_HT-400to600':'ZJets_ht400to600',
        'DYJets_HT-600to800':'ZJets_ht600to800',
        'DYJets_HT-600toInf':'ZJets_ht600toinf',
        'DYJets_HT-800to1200':'ZJets_ht800to1200',
        'DYJets_HT-1200to2500':'ZJets_ht1200to2500',
        'DYJets_HT-2500toInf':'ZJets_ht2500toinf',

        'WJetsToLNu_HT-100to200':'WJets_ht100to200',
        'WJetsToLNu_HT-200to400':'WJets_ht200to400',
        'WJetsToLNu_HT-400to600':'WJets_ht400to600',
        'WJetsToLNu_HT-600to800':'WJets_ht600to800',
        'WJetsToLNu_HT-800to1200':'WJets_ht800to1200',
        'WJetsToLNu_HT-1200to2500':'WJets_ht1200to2500',
        'WJetsToLNu_HT-2500toInf':'WJets_ht2500toinf',

    #'ZtoNuNu'             : ['ZtoNuNu_ht100to200','ZtoNuNu_ht200to400','ZtoNuNu_ht400to600','ZtoNuNu_ht600to800','ZtoNuNu_ht800to1200','ZtoNuNu_ht1200to2500','ZtoNuNu_ht2500toinf'],

        }


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

inbase = '/home/snarayan/home000/ajnlo/'
outbase = inbase+'/merged/'

suffix = ' > /dev/null '
#suffix = ''

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
    f = root.TFile.Open(fpath,'UPDATE')
    t = f.Get('Events')
    n = root.Normalizer();
    n.inWeightName = ''
    n.NormalizeTree(t,t.GetEntries(),xsec)
    f.WriteTObject(t,'Events','Overwrite')
    f.Close()


def merge(shortnames,mergedname):
    for shortname in shortnames:
        t = translate[shortname]
        if t in pds:
            pd = pds[t][0]
            xsec = pds[t][1]
        inpath = inbase+shortname+'_*.root'
        hadd(inpath,'/tmp/%s/split/%s.root'%(user,shortname))
        if xsec>0:
            normalizeFast('/tmp/%s/split/%s.root'%(user,shortname),xsec)
    hadd(['/tmp/%s/split/%s.root'%(user,x) for x in shortnames],'/tmp/%s/merged/%s.root'%(user,mergedname))

d = {
    'a_lo'               : ['GJets_HT-100to200','GJets_HT-200to400','GJets_HT-400to600','GJets_HT-600toInf'],
    'z_lo'               : ['DYJets_HT-100to200',    'DYJets_HT-200to400',    'DYJets_HT-400to600',    'DYJets_HT-600toInf'],
    'w_lo'               : ['WJetsToLNu_HT-100to200',    'WJetsToLNu_HT-200to400',    'WJetsToLNu_HT-400to600',    'WJetsToLNu_HT-600to800',    'WJetsToLNu_HT-800to1200'],
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

