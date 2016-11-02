#!/usr/bin/env python

from sys import argv
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
Load('Drawers','CanvasDrawer')

plot = root.CanvasDrawer()
plot.SetTDRStyle()

masscorr = '' # 'L2L3'

ftemplate = root.TFile('~/public_html/figs/toptagging/datavalidation/v8/templates/tag_top_ecfv8_bdt_pass_hists.root')
hdata = ftemplate.Get('h_fj1MSD%s_Data'%masscorr)
hprong = {}; dhprong = {}; pdfprong = {}; norm = {}
hprong[3] = ftemplate.Get('h_fj1MSD%s_3-prong'%masscorr)
hprong[2] = ftemplate.Get('h_fj1MSD%s_2-prong'%masscorr)
hprong[1] = ftemplate.Get('h_fj1MSD%s_1-prong'%masscorr)

mass = root.RooRealVar("m","m_{SD} [GeV]",40,450)

for iP in xrange(1,4):
  print iP,hprong[iP]
  dhprong[iP] = root.RooDataHist('dh%i'%iP,'dh%i'%iP,root.RooArgList(mass),hprong[iP]) 
  pdfprong[iP] = root.RooHistPdf('pdf%i'%iP,'pdf%i'%iP,root.RooArgSet(mass),dhprong[iP]) 
  norm_ = hprong[iP].Integral()
  norm[iP] = root.RooRealVar('normm%i'%iP,'norm%i'%iP,norm_,0.5*norm_,1.5*norm_)

model = root.RooAddPdf('model','model',root.RooArgList(*[pdfprong[x] for x in [1,2,3]]),root.RooArgList(*[norm[x] for x in [1,2,3]]))

dh_data = root.RooDataHist('dh_data','dh_data',root.RooArgList(mass),hdata)

model.fitTo(dh_data)

frame = mass.frame()
dh_data.plotOn(frame)
hdata.Scale((norm[1].getVal()+norm[2].getVal()+norm[3].getVal())/hdata.Integral())
model.plotOn(frame)
#for c,l in [('ttm',root.RooFit.LineColor(root.kRed))]:
for c,l in [(3,root.RooFit.LineColor(root.kOrange)),(2,root.RooFit.LineColor(6)),(1,root.RooFit.LineColor(8))]:
  model.plotOn(frame,root.RooFit.Components('pdf%i'%c),l)

c = root.TCanvas()
frame.SetTitle('')
frame.Draw()
for h in hprong.values():
  h.SetLineWidth(2)
  h.SetLineStyle(2)
hprong[3].SetLineColor(root.kOrange); hprong[3].Draw('samehist')
hprong[2].SetLineColor(6);            hprong[2].Draw('samehist')
hprong[1].SetLineColor(8);            hprong[1].Draw('samehist')

plot.AddCMSLabel()
plot.SetLumi(12.9); plot.AddLumiLabel(True)
plot.SetCanvas(c)
plot.Draw('~/public_html/figs/toptagging/datavalidation/v8/templates/','simplefit')
