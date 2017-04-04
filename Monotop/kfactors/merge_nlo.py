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


VERBOSE=False

user = environ['USER']
system('mkdir -p /tmp/%s/split'%user) # tmp dir
system('mkdir -p /tmp/%s/merged'%user) # tmp dir

inbase = '/home/snarayan/home000/ajnlo/'
outbase = inbase+'/merged/'

N_EVENTS_FILE = 10000

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

def normalizeFast(fpath,nfiles):
    if 'unlops' not in fpath:
        nevts = nfiles * N_EVENTS_FILE
    else:
        nevts = nfiles
    PInfo(sname,'normalizing %s (%s) ...'%(fpath,nfiles))
    f = root.TFile.Open(fpath,'UPDATE')
    t = f.Get('Events')
    n = root.Normalizer();
    n.inWeightName = 'evtweight'
    n.NormalizeTree(t,nevts,1.)
    f.WriteTObject(t,'Events','Overwrite')
    f.Close()


def merge(shortnames,mergedname):
    for shortname in shortnames:
        inpath = inbase+shortname+'_*.root'
        nfiles = len(glob(inpath))
        if 'unlops' in shortname:
            nfiles /= 1.e9
        hadd(inpath,'/tmp/%s/split/%s.root'%(user,shortname))
        normalizeFast('/tmp/%s/split/%s.root'%(user,shortname),nfiles)
    hadd(['/tmp/%s/split/%s.root'%(user,x) for x in shortnames],'/tmp/%s/merged/%s.root'%(user,mergedname))

d = {
        'a_nlo' : ['unlops'],
        'z_nlo' : ['z12j'],
        'w_nlo' : ['w12j'],
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

