#!/usr/bin/env python

from os import getenv
from PandaCore.Tools.ConfigBuilder import *
import subprocess
import sys
import argparse
from glob import glob
from re import sub

outdir = getenv('SUBMIT_OUTDIR')
parser = argparse.ArgumentParser(description='check missing files')
parser.add_argument('--infile',type=str)
parser.add_argument('--outfile',type=str)
parser.add_argument('--outdir',type=str,default=outdir)
args = parser.parse_args()
outdir = args.outdir

class Output:
  def __init__(self,name):
    self.name = name
    self.total = 0
    self.done = 0
  def add(self,isComplete):
    self.total += 1
    if isComplete:
      self.done += 1
  def __str__(self):
    if self.total==0:
      return ''
    s = '%-50s'%self.name
    frac = 1.*self.done/self.total
    s += '\t[\033[0;42m'
    switched = False
    for i in xrange(100):
      if i/100.<frac or switched:
        s += ' '
      else:
        s += '\033[0;41m '
        switched = True
    s += '\033[0m] %i/%i (%.2f%%)'%(self.done,self.total,100*frac)
    return s

sys.argv=[]

processedfiles = []

locks = glob(outdir+'/locks/*lock')
for lock in locks:
	flock = open(lock)
	for l in flock:
		processedfiles.append(l.strip())

outputs = {}
data = Output('Data')
mc = Output('MC')

all_samples = read_sample_config(args.infile)
filtered_samples = {}
outfile = open(args.outfile,'w')
for name in sorted(all_samples):
	sample = all_samples[name]
	out_sample = DataSample(name,sample.dtype,sample.xsec)

	base_name = sub('_[0-9]+$','',name)
	if base_name not in outputs:
		outputs[base_name] = Output(base_name)
	output = outputs[base_name]

	for f in sample.files:
		found = (f in processedfiles)
		if not found:
			out_sample.add_file(f)
		output.add(found)
		if sample.dtype=='MC':
			mc.add(found)
		else:
			data.add(found)

	if len(out_sample.files)>0:
		filtered_samples[name] = out_sample

keys = sorted(filtered_samples)
for k in keys:
	sample = filtered_samples[k]
	configs = sample.get_config(-1)
	for c in configs:
		outfile.write(c)

for n in sorted(outputs):
  print str(outputs[n])
print str(data)
print str(mc)

outfile.close()

