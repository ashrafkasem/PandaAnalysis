#!/usr/bin/env python

from glob import glob
from os import stat,getenv
from PandaCore.Tools.process import *
from re import sub
import sys
import argparse

print 'starting'

parser = argparse.ArgumentParser(description='make config file')
parser.add_argument('--outfile',type=str)
parser.add_argument('--include',nargs='+',type=str,default=None)
parser.add_argument('--exclude',nargs='+',type=str,default=None)
args = parser.parse_args()

proddirs = getenv('PRIVATE_PRODDIR').split(':')

listOfFiles=[]
for d in proddirs:
  print 'searching in',d+'/*/panda_*.root'
  listOfFiles += glob(d+'/*/panda_*.root')

def checkDS(nickname,include,exclude):
  included=False
  if include:
    for i in include:
      if i in nickname:
        included=True
        break
  else:
    included=True
  excluded=False
  if exclude:
    for e in exclude:
      if e in nickname:
        excluded=True
        break
  else:
    excluded=False
  return (included and not(excluded))

print 'found %i files'%len(listOfFiles)
print 'writing to',args.outfile
cfgFile = open(args.outfile,'w')

couldNotFind = []

file_counter=0
for f in sorted(listOfFiles):
  if stat(f).st_size==0:
    continue
  ff = f.split('/')
  pd = ff[-2]
  nickname = pd
  if checkDS(nickname,args.include,args.exclude):
    nickname += '_%i'%(file_counter)
    cfgFile.write('{0:<60} {2:<6} {3:<15} {1:<150}\n'.format(nickname,f,'MC',1)) 
    file_counter += 1

if len(couldNotFind)>0: 
  print 'could not find:'
for pd in couldNotFind:
  print '\t',pd
