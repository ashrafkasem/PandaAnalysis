#!/usr/bin/env python

from sys import argv,exit
from os import getenv
from array import array
from math import sqrt

argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
Load('Drawers','HistogramDrawer')

root.gROOT.LoadMacro('triggerFunc.C')
fitFunc = root.TF1('fitFunc',root.ErfCB,60,900,5)
initParams = [80,1,0.064,1.72,1]
for iP in xrange(5):
  fitFunc.SetParameter(iP,initParams[iP])
  #fitFunc.SetParameter(iP,2)

plot = root.HistogramDrawer()
plot.SetLumi(12.9)
plot.SetTDRStyle()

counter=0

fIn = root.TFile('/home/snarayan/home000/panda/vbf_v0/SingleElectron.root')
events = fIn.Get('events')
hOne = root.TH1F('hone','one',1,0.,2.)

baseCut = 'metFilter==1 && ((trigger&2)!=0) && nLooseElectron==1 && looseLep1IsTight && nLooseMuon==0 && nJet>1 && jet1Pt>%f && jet2Pt>%f && jet12DEta>3'
monobaseCut = 'metFilter==1 && ((trigger&2)!=0) && nLooseElectron==1 && looseLep1IsTight && nLooseMuon==0 && nJet>0 && jet1Pt>%f && fabs(jet1Eta)<2.5'
#monobaseCut = '((trigger&2)!=0) && nLooseElectron>0'
metCut  = '(trigger&1)!=0'

def getHist(tree,cut,bins):
  global counter
  N = len(bins)-1
  h = root.TH1F('h%i'%counter,'h%i'%counter,N,bins)
  tree.Draw('metnomu>>h%i'%counter,cut)
  counter += 1
  return h

def getMETHist(tree,cut,bins):
  global counter
  N = len(bins)-1
  h = root.TH1F('h%i'%counter,'h%i'%counter,N,bins)
  tree.Draw('pfmet>>h%i'%counter,cut)
  counter += 1
  return h

def runEff(label,fn=getHist,doFit=False,pt1=100,pt2=40):
  bins = array('f',[60,80,100,120,140,160,180,200,220,240,270,300,400,900])
  baseCut_ = baseCut%(pt1,pt2)
  if pt2==0:
    #baseCut_ = monobaseCut
    baseCut_ = monobaseCut%(pt1)
  hBase = fn(events,baseCut_,bins)
  hMET = fn(events,tAND(baseCut_,metCut),bins)
  x = []; central = []; up = []; down = []; zero = []
  for iB in xrange(1,len(bins)):
    passed=hMET.GetBinContent(iB)
    total=hBase.GetBinContent(iB)
    try:
      eff = passed/total
    except ZeroDivisionError:
      eff = 0
    try:
#      err = eff*sqrt(1/total+1/passed)
#      err = sqrt(pow(passed,2)/pow(total,3)+passed/pow(total,2))
#      err = sqrt(passed*(1-passed/total))/total 
       errUp = root.TEfficiency.ClopperPearson(int(total),int(passed),0.68,True)
       errDown = root.TEfficiency.ClopperPearson(int(total),int(passed),0.68,False)
    except ZeroDivisionError:
       err = eff
       errUp = eff
       errDown = eff
    hMET.SetBinContent(iB,eff)
    hMET.SetBinError(iB,(errUp-errDown)/2)
    x.append(hMET.GetBinCenter(iB))
    central.append(eff)
    up.append(errUp-eff)
    down.append(eff-errDown)
#    up.append(min(1-eff,err))
#    down.append(min(eff,err))
    zero.append(0)
  N = len(x)
  x = array('f',x); central = array('f',central); zero=array('f',zero)
  up=array('f',up); down=array('f',down)
  errs = root.TGraphAsymmErrors(N,x,central,zero,zero,down,up)
  errs.SetLineWidth(2)
  if fn==getHist:
    hMET.GetXaxis().SetTitle('U')
  else:
    hMET.GetXaxis().SetTitle('MET')
  hMET.GetYaxis().SetTitle('Efficiency')
  hMET.SetMaximum(1)
  plot.AddCMSLabel()
  plot.AddLumiLabel(True)
  if doFit:
    hMET.Fit(fitFunc)
    plot.AddPlotLabel('#mu=%.3f, #sigma=%.3f'%(fitFunc.GetParameter(0),fitFunc.GetParameter(1)),0.5,0.5,False,42,0.05,11)
  if pt2:
    plot.AddPlotLabel('p_{T}^{1}>%.0f, p_{T}^{2}>%.0f'%(pt1,pt2),0.5,0.4,False,42,0.05,11)
  else:
    plot.AddPlotLabel('p_{T}^{1}>%.0f'%(pt1),0.5,0.4,False,42,0.05,11)
  for iB in xrange(1,len(bins)):
    hMET.SetBinError(iB,0.00001)
  fitFunc.SetLineColor(root.kBlue)
  plot.AddHistogram(hMET,'Data',root.kData)
  plot.AddAdditional(errs,'p')
  if doFit:
    plot.AddAdditional(fitFunc,'l')
  plot.Draw('~/public_html/figs/vbf/trigger/',label)
  plot.Reset()

fit=False

runEff('vbf_baseline',getMETHist,fit) 
runEff('vbf_pt1_40_pt2_40',getMETHist,fit,40,40) 
runEff('vbf_pt1_20_pt2_20',getMETHist,fit,20,20) 

#runEff('mono_baseline',getMETHist,fit,100,0) 
#runEff('mono_pt1_40',getMETHist,fit,40,0) 
#runEff('mono_pt1_20',getMETHist,fit,20,0) 

