#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv,exit
import argparse

basedir = getenv('PANDA_FLATDIR')
outdir = basedir+'/sf/'

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--indir',metavar='indir',type=str,default=basedir)
parser.add_argument('--outdir',metavar='outdir',type=str,default=outdir)
parser.add_argument('--sel',metavar='sel',type=str,default='tag')
parser.add_argument('--cat',metavar='cat',type=str)
args = parser.parse_args()

outdir = args.outdir
basedir = args.indir
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
tcut = root.TCut
Load('PandaAnalysisFlat','SFTreeBuilder')

### SET GLOBAL VARIABLES ###
lumi = 12918.
logy=False
nlo = 'sf_ewkV*sf_qcdV'
if args.sel=='mistag':
  cut = 'nFatjet==1 && fj1Pt>250 && fj1MaxCSV<0.46 && nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && UWmag>250 && isojetNBtags==0'
  weight = '%f*normalizedWeight*sf_pu*sf_lep*%s*sf_sjbtag0*sf_btag0*sf_tt*sf_metTrig'%(lumi,nlo)
  outfileprefix = 'mistag_'
elif args.sel=='photon':
  cut = 'nFatjet==1 && fj1Pt>250 && nLooseLep==0 && nLoosePhoton==1 && loosePho1IsTight==1 && nTau==0 && UAmag>250'
  weight = '%f*normalizedWeight*sf_pu*sf_lep*%s*sf_tt*sf_phoTrig*0.93'%(lumi,nlo)
  outfileprefix = 'photon_'
else:
  cut = 'nFatjet==1 && fj1Pt>250 && fj1MaxCSV>0.46 && nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && UWmag>250 && isojetNBtags==1'
  weight = '%f*normalizedWeight*sf_pu*sf_lep*%s*sf_sjbtag1*sf_btag1*sf_tt*sf_metTrig'%(lumi,nlo)
  outfileprefix = 'tag_'

if args.cat=='3':
  cut = tAND(cut,'fj1IsMatched==1&&fj1GenSize<1.44')
  samples = ['TTbar','SingleTop']
  outname='prong3'
elif args.cat=='2':
  cut = tAND(cut,'fj1IsMatched==0||fj1GenSize>1.44')
  samples = ['TTbar','Diboson','SingleTop']
  outname='prong2'
elif args.cat=='1':
  samples = ['WJets','QCD']
  outname='prong1'
elif args.cat=='gjets':
  samples = ['GJets','QCD']
  outname='gjets'
elif args.cat=='data':
  if args.sel == 'photon':
    cut = tAND(cut,'(trigger&4)!=0')
    samples = ['SinglePhoton']
  else:
    cut = tAND(cut,'(trigger&1)!=0')
    samples = ['MET']
  weight = '1'
  outname='data'
else:
  PError(sname,'Category %s not recognized'%args.cat)
  exit(1)

# load the builder
skimmer = root.SFTreeBuilder()

skimmer.AddTagger('fj1Tau32','tau32')
skimmer.AddTagger('fj1Tau32SD','tau32SD')
skimmer.AddTagger('top_ecfv6_bdt','top_ecfv6_bdt')
skimmer.AddTagger('top_ecfv8_bdt','top_ecfv8_bdt')
skimmer.cutFormula = cut
skimmer.weightFormula = weight
for s in samples:
  skimmer.AddInputFile(basedir+'/'+s+'.root')
skimmer.SetOutputFile(outdir+'/%s_%s.root'%(outfileprefix,outname),outname)
skimmer.Run()
skimmer.Terminate()

