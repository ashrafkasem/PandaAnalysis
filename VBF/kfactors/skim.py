#!/usr/bin/env python

from re import sub
from sys import argv,exit
import argparse
argv=[]

from os import system,getenv
from PandaCore.Tools.Load import *
import ROOT as root

outdir = '~/home000/store/kfactors/'
Load('PandaAnalysisFlat','GenAnalyzer')

def fn(infiles,outname,order,ptype):

  skimmer = root.GenAnalyzer()
 
  skimmer.order = order
  skimmer.processType = ptype


  skimmer.SetOutputPath(outdir+'skimmed',outname)
  files = []
  for f in infiles:
    fin = root.TFile(outdir+f[0])
    files.append(fin)
    if order==root.GenAnalyzer.kNLO:
      tree = fin.FindObjectAny("genEvents")
    else:
      tree = fin.FindObjectAny("events")
    skimmer.AddInput(tree,f[1],f[2])

  skimmer.Run()
  skimmer.Terminate()


#fn([('lo_wlo.root',1.5,'lo_wlo')],'w_lo.root',root.GenAnalyzer.kLO,root.GenAnalyzer.kW)
#fn([('lo_zlo.root',1.,'lo_zlo')],'z_lo.root',root.GenAnalyzer.kLO,root.GenAnalyzer.kZ)
fn([('nn012j_5f_pt150_NLO_FXFX_genkine.root',16.49*3,'nlo_z1'),
    ('nn012j_5f_pt80_NLO_FXFX_genkine.root',90.81*3,'nlo_z0')],
    'z_nlo.root',
    root.GenAnalyzer.kNLO,
    root.GenAnalyzer.kZ)
'''
fn([('nl012j_5f_pt80_NLO_FXFX_genkine.root',388.3*3,'nlo_w0'),
    ('nl012j_5f_pt150_NLO_FXFX_genkine.root',65.8*3,'nlo_w1')],
    'w_nlo.root',
    root.GenAnalyzer.kNLO,
    root.GenAnalyzer.kW)
'''
