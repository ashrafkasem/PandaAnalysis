#!/usr/bin/env python

from os import getenv
from PandaCore.Tools.process import *
import subprocess
import sys
import argparse

outdir = getenv('SUBMIT_OUTDIR')
parser = argparse.ArgumentParser(description='check missing files')
parser.add_argument('--infile',type=str,default=None)
parser.add_argument('--outfile',type=str,default=None)
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
    #s = self.name+'\n'
    frac = 1.*self.done/self.total
    s += '\t[\033[0;42m'
    #s += '\t[\033[0;32m'
    switched = False
    for i in xrange(100):
      if i/100.<frac or switched:
        s += ' '
      else:
        s += '\033[0;41m '
        switched = True
    s += '\033[0m] %i/%i (%.2f%%)'%(self.done,self.total,100*frac)
    #s += '\033[0m] %i/%i (%.2f%%)'%(self.done,self.total,100*frac)
    return s

sys.argv=[]

foundfiles = []
cmd = 'ls %s'%(outdir)
print cmd
for line in subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout.readlines():
  foundfiles.append(line.strip())

outputs = {}
data = Output('Data')
mc = Output('MC')

infile = open(args.infile)
nmissing=0; ndata=0
if args.outfile:
  outfile = open(args.outfile,'w')
else:
  outfile = open(args.infile.replace('.cfg','_missing.cfg'),'w')
for line in infile:
  shortname = line.split()[0] + '.root'
  samplename = '_'.join(shortname.split('_')[:-1])
  if samplename not in outputs:
    output = Output(samplename)
    outputs[samplename] = output
  else:
    output = outputs[samplename]
  found=False
  for f in foundfiles:
    if shortname==f:
      found=True
      break
  if not found:
    outfile.write(line)
    if 'Data' in line:
      data.add(False)
    else:
      mc.add(False)
    output.add(False)
  else:
    if 'Data' in line:
      data.add(True)
    else:
      mc.add(True)
    output.add(True)

for n in sorted(outputs):
  print str(outputs[n])
print str(data)
print str(mc)

outfile.close()
infile.close()

