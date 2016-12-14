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

ptBins = [100,120,160,200,250,300,350,400,500,600,1000]
nPtBins = len(ptBins)-1
ptBins = array('f',ptBins)

plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.AddCMSLabel()
plot.Logy(True)
plot.SetAbsMin(0.0001)
plot.SetTDRStyle()
plot.InitLegend()

plotr = root.HistogramDrawer()
plotr.SetLumi(12.9)
plotr.SetRatioStyle()
plotr.AddCMSLabel()
plotr.AddLumiLabel()
plotr.InitLegend(.15,.6,.5,.8)

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

def plotDist(V,dists):
  if V=='Z':
    tlo = tzlo
    tnlo = tznlo
  else:
    tlo = twlo
    tnlo = twnlo
  toreturn = []
  for d in dists:
    hlo = getDist(tlo,d[0],d[1],d[2])
    hnlo = getDist(tnlo,d[0],d[1],d[2])
    toreturn.append((hlo,hnlo))
    plot.AddHistogram(hlo,'%s LO'%(V),root.kSignal2)
    plot.AddHistogram(hnlo,'%s NLO'%(V),root.kExtra2)
    if len(d)<4 or d[3]==None:
      plot.Draw(args.outdir,V+'_'+d[0])
    else:
      plot.Draw(args.outdir,V+'_'+d[3])
    plot.Reset()
  return toreturn

'''
 plotr.AddHistogram(hsigewk,'Signal region',root.kExtra2,root.kBlack)
  if V=='W':
    plotr.AddHistogram(hmuewk,'Single lep CRs',root.kExtra1)
  #  plotr.AddHistogram(heleewk,'Single e CR',root.kExtra4)
  else:
    plotr.AddHistogram(hmuewk,'Dilep CRs',root.kExtra1)
  #  plotr.AddHistogram(heleewk,'Dielectron CR',root.kExtra4)
  plotr.AddAdditional(hsigerr,'e2')
  plotr.AddAdditional(hmuerr,'e2')
  plotr.Draw(args.outdir,label)
  plotr.Reset()

def plotWZ(label):
  totalcut = tTIMES(tAND(vbfsel,cuts['signal']),weights['signal'])
  hqcdz = getMETHist(tqcdzvv,totalcut)
  hqcdw = getMETHist(tqcdw,totalcut)
  hewkz = getMETHist(tewkzvv,totalcut)
  hewkw = getMETHist(tewkw,totalcut)
  hewkw.Divide(hqcdw)
  hewkz.Divide(hqcdz)
  for h in [hewkw,hewkz]:
    h.GetYaxis().SetTitle('EWK V/QCD V')
  hwerr = hewkw.Clone()
  hwerr.SetFillStyle(3004)
  hwerr.SetFillColorAlpha(root.kGreen+1,0.5)
  hwerr.SetLineWidth(0)
  hzerr = hewkz.Clone()
  hzerr.SetFillStyle(3005)
  hzerr.SetFillColorAlpha(root.kCyan+1,0.5)
  hzerr.SetLineWidth(0)
  plotr.AddHistogram(hewkz,'EWK Z/QCD Z',root.kExtra2)
  plotr.AddHistogram(hewkw,'EWK W/QCD W',root.kExtra3)
  plotr.AddAdditional(hzerr,'e2')
  plotr.AddAdditional(hwerr,'e2')
  plotr.Draw(args.outdir,label)
  plotr.Reset()

'''

plotDist('Z',[('vpt',recoilBins,'p_{T}^{V} [GeV]')])
plotDist('W',[('vpt',recoilBins,'p_{T}^{V} [GeV]')])

