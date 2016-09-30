#!/usr/bin/env python

import cPickle as pickle
from sys import argv

nFiles = len(argv)-1

badEventFilter=set([])
for i in xrange(1,nFiles):
  fTxt = open(argv[i+1])
  lines = list(fTxt)
  print argv[i+1],len(lines)
  for line in lines:
    ll = line.strip().split(':')
    badEventFilter.add((int(ll[0]),int(ll[1]),int(ll[2])))

with open(argv[1],'wb') as pklFile:
  pickle.dump(badEventFilter,pklFile,-1)
