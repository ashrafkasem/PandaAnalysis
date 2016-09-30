#!/usr/bin/env python

import cPickle as pickle
from sys import argv
from re import sub

fIn = open(argv[1],'rb')
allFilter = pickle.load(fIn)
nEvts = len(allFilter)
nPerFile = int(nEvts/10)
evtCounter=0
fileCounter=0
subFilter = set([])
for e in allFilter:
  evtCounter += 1
  if (evtCounter%nPerFile)==0 and fileCounter<9:
    print fileCounter,len(subFilter),e
    with open(sub('.pkl','_%i.pkl'%fileCounter,argv[1]),'wb') as fOut:
      pickle.dump(subFilter,fOut,-1)
    del subFilter
    subFilter = set([])
    fileCounter += 1
  subFilter.add(e)

with open(sub('.pkl','_9.pkl',argv[1]),'wb') as fOut:
  pickle.dump(subFilter,fOut,-1)
