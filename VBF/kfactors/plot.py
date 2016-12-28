#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = '/home/snarayan/home000/store/kfactors/skimmed/'
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
args = parser.parse_args()
sname=argv[0]
argv=[]

import ROOT as root
from ROOT import gROOT
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
from array import array
from math import sqrt
Load('Drawers','HistogramDrawer')

### DEFINE REGIONS ###

recoilBins = [200,250,300,350,400,500,600,1000]
nRecoilBins = len(recoilBins)-1
recoilBins = array('f',recoilBins)

ptBins = [100,120,160,200,250,300,350,400,450,500,550,600,650,700,800,900,1000,1200]
nPtBins = len(ptBins)-1
ptBins = array('f',ptBins)

plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.AddCMSLabel()
plot.Logy(True)
#plot.SetAbsMin(0.0001)
plot.InitLegend()

plotr = root.HistogramDrawer()
plotr.SetRatioStyle()
plotr.AddCMSLabel()
#plotr.InitLegend(.15,.6,.5,.8)
plotr.InitLegend()

counter=0

fzlo = root.TFile(baseDir+'z_lo.root'); tzlo = fzlo.Get('events')
fznlo = root.TFile(baseDir+'z_nlo.root'); tznlo = fznlo.Get('events')
fwlo = root.TFile(baseDir+'w_lo.root'); twlo = fwlo.Get('events')
fwnlo = root.TFile(baseDir+'w_nlo.root'); twnlo = fwnlo.Get('events')

ctmp = root.TCanvas()

def getDist(tree,var,bins,xlabel,cut='1==1'):
  global counter
  ctmp.cd()
  if len(bins)==3:
    h = root.TH1D('h%i'%counter,'h%i'%counter,bins[0],bins[1],bins[2])
    scale=False
  else:
    h = root.TH1D('h%i'%counter,'h%i'%counter,len(bins)-1,bins)
    scale=True
  h.GetXaxis().SetTitle(xlabel)
  h.GetYaxis().SetTitle('')
  tree.Draw('%s>>h%i'%(var,counter),'weight*(%s)'%(cut))
  if scale:
    h.Scale(1,'width')
  counter += 1
  h.SetFillStyle(0)
  return h

def plotDist(V,dists,cut):
  if V=='Z':
    tlo = tzlo
    tnlo = tznlo
  else:
    tlo = twlo
    tnlo = twnlo
  toreturn = []
  for d in dists:
    hlo = getDist(tlo,d[0],d[1],d[2],cut)
    hnlo = getDist(tnlo,d[0],d[1],d[2],cut)
    toreturn.append((hlo,hnlo))
    plot.AddHistogram(hlo,'%s LO'%(V),root.kSignal2)
    plot.AddHistogram(hnlo,'%s NLO'%(V),root.kExtra2)
    if len(d)<4 or d[3]==None:
      plot.Draw(args.outdir,V+'_'+d[0])
    else:
      plot.Draw(args.outdir,V+'_'+d[3])
    plot.Reset()
    plot.AddCMSLabel()
  return toreturn

def plotKFactors(V,hists,name):
  # hists is a list of tuples (hlo, hnlo, label)
  counter=0
  for hlo,hnlo,label in hists:
    hratio = hnlo.Clone()
    hratio.Divide(hlo)
    if counter==0:
      hratio.SetMaximum(2); hratio.SetMinimum(0)
    plotr.AddHistogram(hratio,label,root.kExtra1+counter)
    hratioerr = hratio.Clone()
    hratioerr.SetFillStyle(3004)
    hratioerr.SetFillColorAlpha(root.kBlack,0.5)
    hratioerr.SetLineWidth(0)
    plotr.AddAdditional(hratioerr,'e2')
    counter += 1
  plotr.Draw(args.outdir,V+'_'+name)
  plotr.Reset()
  plotr.AddCMSLabel()

hmono = plotDist('Z',[('vpt',ptBins,'p_{T}^{V} [GeV]','vpt_monojet')],'njet>0 && jet1pt>100')[0]
hdi   = plotDist('Z',[('vpt',ptBins,'p_{T}^{V} [GeV]','vpt_dijet')],'njet>1 && jet1pt>80 && jet2pt>40 && jet1eta*jet2eta<0')[0]
hvbf  = plotDist('Z',[('vpt',ptBins,'p_{T}^{V} [GeV]','vpt_vbf')],'njet>1 && jet1pt>80 && jet2pt>40 && jet1eta*jet2eta<0 && mjj>1100')[0]
plotKFactors('Z',[(hmono[0],hmono[1],'Monojet'),
                  (hdi[0],hdi[1],'Dijet'),
                  (hvbf[0],hvbf[1],'VBF')],'kfactor_ptV')
#plotDist('W',[('vpt',recoilBins,'p_{T}^{V} [GeV]')])

