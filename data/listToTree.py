#!/usr/bin/env python

from sys import argv
from ROOT import TFile,TTree
from array import array

nFiles = len(argv)-1

badEventFilter=set([])
for i in xrange(1,nFiles):
  fTxt = open(argv[i+1])
  lines = list(fTxt)
  print argv[i+1],len(lines)
  for line in lines:
    ll = line.strip().split(':')
    badEventFilter.add((int(ll[0]),int(ll[1]),int(ll[2])))

fOut = TFile(argv[1],'RECREATE')
events = TTree('events','events')
r = array('I',[0])
l = array('I',[0])
e = array('I',[0])
events.Branch('run',r,'run/i')
events.Branch('lumi',l,'lumi/i')
events.Branch('evt',e,'evt/i')

for ev in badEventFilter:
  try:
    r[0] = ev[0]
    l[0] = ev[1]
    e[0] = ev[2]
    events.Fill()
  except:
    print ev

fOut.WriteTObject(events)
