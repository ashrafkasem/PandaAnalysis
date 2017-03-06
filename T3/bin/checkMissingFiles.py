#!/usr/bin/env python

from os import getenv
from PandaCore.Tools.job_management import *
import subprocess
import sys
import argparse
from glob import glob
from re import sub
import cPickle as pickle
from itertools import chain 

outdir = getenv('SUBMIT_OUTDIR')
workdir = getenv('SUBMIT_WORKDIR')
parser = argparse.ArgumentParser(description='check missing files')
parser.add_argument('--infile',type=str,default=None)
parser.add_argument('--outfile',type=str,default=None)
parser.add_argument('--outdir',type=str,default=outdir)
parser.add_argument('--force',action='store_true')
parser.add_argument('--nfiles',type=int,default=-1)
args = parser.parse_args()
outdir = args.outdir

if not args.infile:
  args.infile = workdir+'/local_all.cfg'
if not args.outfile:
  args.outfile = workdir+'/local.cfg'


class Output:
  def __init__(self,name):
    self.name = name
    self.total = 0
    self.done = 0
    self.idle = 0
    self.running = 0
    self.missing = 0
  def add(self,state):
    self.total += 1
    if state=='done':
        self.done += 1
    elif state=='running':
        self.running += 1
    elif state=='idle':
        self.idle += 1
    elif state=='missing':
        self.missing += 1
  def __str__(self):
    WIDTH=50
    if self.total==0:
      return ''
    s = '%-80s'%self.name
    d_frac = 1.*WIDTH*self.done/self.total
    r_frac = 1.*WIDTH*(self.done+self.running)/self.total
    i_frac = 1.*WIDTH*(self.idle+self.done+self.running)/self.total
    s += '\t[\033[0;42m'
    state = 0
    for i in xrange(WIDTH):
        if i>=d_frac:
            s += '\033[0;44m'
        if i>=r_frac:
            s += '\033[0;47m'
        if i>=i_frac:
            s += '\033[0;41m'
        s += ' '
    s += '\033[0m] %i/%i (%.2f%%)'%(self.done,self.total,d_frac)
    return s


# determine what files have been processed and logged as such
processedfiles = []
locks = glob(outdir+'/locks/*lock')
for lock in locks:
    flock = open(lock)
    for l in flock:
        processedfiles.append(l.strip())


# determine what samples from previous resubmissions are still running
running_samples = []
idle_samples = []
with open(workdir+'/submission.pkl','rb') as fpkl:
  submissions = pickle.load(fpkl)
for s in submissions:
  results = s.query_status()
  running_samples += results['running']
  idle_samples += results['idle']

running_files = list(chain.from_iterable([x.files for x in running_samples]))
idle_files = list(chain.from_iterable([x.files for x in idle_samples]))


# for fancy display
outputs = {}
data = Output('Data')
mc = Output('MC')


all_samples = read_sample_config(args.infile)
filtered_samples = {}
merged_samples = {}
outfile = open(args.outfile,'w')
for name in sorted(all_samples):
    sample = all_samples[name]
    out_sample = DataSample(name,sample.dtype,sample.xsec)

    base_name = sub('_[0-9]+$','',name)
    if base_name not in outputs:
        outputs[base_name] = Output(base_name)
    output = outputs[base_name]
    if base_name not in merged_samples:
        merged_samples[base_name] = DataSample(base_name,sample.dtype,sample.xsec)
    merged_sample = merged_samples[base_name]

    to_resubmit = []

    for f in sample.files:
        state = 'missing'
        if f in processedfiles:
            state = 'done'
        elif f in running_files:
            state = 'running'
        elif f in idle_files:
            state = 'idle'

        if state=='missing' or (args.force and state!='done'):
            # if '750' in f:
            #     print '|%s|'%f
            out_sample.add_file(f)
            merged_sample.add_file(f)

        output.add(state)
        if sample.dtype=='MC':
            mc.add(state)
        else:
            data.add(state)

    if len(out_sample.files)>0:
        filtered_samples[name] = out_sample

if args.nfiles<0:
    keys = sorted(filtered_samples)
    for k in keys:
        sample = filtered_samples[k]
        if len(sample.files)==0:
            continue
        configs = sample.get_config(-1)
        for c in configs:
            outfile.write(c)
else:
    keys = sorted(merged_samples)
    counter=0
    for k in keys:
        sample = merged_samples[k]
        if len(sample.files)==0:
            continue
        configs = sample.get_config(args.nfiles,suffix='_%i')
        for c in configs:
            outfile.write(c%(counter,counter))
            counter += 1

for n in sorted(outputs):
  print str(outputs[n])
print str(data)
print str(mc)

outfile.close()

