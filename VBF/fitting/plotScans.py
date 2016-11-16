#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
from glob import glob
import argparse
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--var1',type=str)
parser.add_argument('--var2',type=str)
args = parser.parse_args()
figsdir = args.outdir
var1=args.var1; var2=args.var2
basedir = getenv('PANDA_VBFSCAN')
argv=[]
labels = {
    'eta':'jjDEta',
    'phi':'jjDPhi',
    'mjj':'mjj',
    'jet1':'jot1Pt',
    'jet2':'jot2Pt',
    'dphi':'fabsminJetMetDPhiwithendcap',
    }
titles = {
    'eta':'min #Delta#eta(j_{1},j_{2})',
    'phi':'max #Delta#phi(j_{1},j_{2})',
    'mjj':'min M(j_{1},j_{2}) [GeV]',
    'jet1':'min p_{T}(j_{1}) [GeV]',
    'jet2':'min p_{T}(j_{2}) [GeV]',
    'dphi':'min #Delta#phi(j,MET)',
    }

import ROOT as root
from PandaCore.Tools.Load import *

##Color palette
'''
ncontours = 999;
stops = [     0.10,     0.2000,      0.40,   0.600,     0.800,    1.0000]
red   = [255./255.,  169./255.,  91./255., 36./255.,   8./255.,   6./255.]
green = [255./255.,  194./255., 149./255., 99./255.,  46./255.,   0./255.]
blue  = [255./255.,  243./255., 240./255.,240./255., 240./255., 233./255.]
stopsArray = array('d',stops)
redArray   = array('d',red)
greenArray = array('d',green)
blueArray  = array('d',blue)
root.TColor.CreateGradientColorTable(len(stops), stopsArray, redArray, greenArray, blueArray, ncontours);
'''
root.gStyle.SetNumberContours(999);
root.gStyle.SetPalette(root.kBird)
root.gStyle.SetPaintTextFormat(".2g")


Load('Drawers','CanvasDrawer')

plot = root.CanvasDrawer()
plot.SetTDRStyle()
root.gStyle.SetPadRightMargin(0.15)

c = root.TCanvas()

deta = set([]); dphi = set([])
listoffiles = glob(basedir+'/higgsCombine%s*%s*.root'%(labels[var1],labels[var2]))
for f in listoffiles:
  try:
    fname = f.split('/')[-1]
    fname = fname.replace('higgsCombine','').replace('.Asymptotic.mH120.root','')
    etacut = float(fname.split('_')[1])
    phicut = float(fname.split('_')[-1])
    deta.add(etacut)
    dphi.add(phicut)
  except IndexError:
    pass

def xformSet(s):
  s_ = sorted(list(s))
  l = [(3*s_[0]-s_[1])/2]
  l += [(s_[x]+s_[x+1])/2 for x in xrange(len(s_)-1)]
  l += [(3*s_[-1]-s_[-2])/2]
  return l

etaarray = array('f',xformSet(deta))
phiarray = array('f',xformSet(dphi))
h2 = root.TH2F('h2','h2',len(etaarray)-1,etaarray,len(phiarray)-1,phiarray)
h2x = h2.GetXaxis(); h2y = h2.GetYaxis();
h2x.SetTitle(titles[var1])
h2y.SetTitle(titles[var2])
h2.GetZaxis().SetTitle('Expected limit')

for e in deta:
  for p in dphi:
    try:
      if var1=='mjj':
        fname = basedir+'/higgsCombine%s_%i__%s_%.2f.Asymptotic.mH120.root'%(labels[var1],int(e),labels[var2],p)
      elif var2=='mjj':
        fname = basedir+'/higgsCombine%s_%.2f__%s_%i.Asymptotic.mH120.root'%(labels[var1],e,labels[var2],int(p))
      else:
        fname = basedir+'/higgsCombine%s_%.2f__%s_%.2f.Asymptotic.mH120.root'%(labels[var1],e,labels[var2],p)
      f = root.TFile(fname)
      t = f.Get('limit')
      t.GetEntry(2)
      h2.SetBinContent(h2x.FindBin(e),h2y.FindBin(p),t.limit)
    except:
      pass

c.Clear(); c.cd()
h2.Draw('colz text')
plot.AddCMSLabel(.16,.94)
plot.SetLumi(12.9); plot.AddLumiLabel(True)
plot.SetCanvas(c)
plot.Draw(args.outdir+'/','optimized_scan_'+var1+var2)

