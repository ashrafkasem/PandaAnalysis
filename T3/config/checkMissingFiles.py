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

sys.argv=[]

foundfiles = []
cmd = 'ls %s'%(outdir)
print cmd
for line in subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout.readlines():
  foundfiles.append(line.strip())

infile = open(args.infile)
nmissing=0; ndata=0
if args.outfile:
  outfile = open(args.outfile,'w')
else:
  outfile = open(args.infile.replace('.cfg','_missing.cfg'),'w')
for line in infile:
  shortname = line.split()[0] + '.root'
  found=False
  for f in foundfiles:
    if shortname==f:
      found=True
      break
  if not found:
    nmissing+=1
    if 'Data' in line:
      ndata+=1
    outfile.write(line)
print 'Missing:',nmissing,'of which',ndata,'is data'
outfile.close()
infile.close()

