#!/usr/bin/env python

from sys import argv
argv=[]

import ROOT as root

ftemplate = root.TFile('~/public_html/figs/toptagging/datavalidation/v7/templates/fail_hists.root')
hdata = ftemplate.Get('h_fj1MSD_Data')
hmatched = ftemplate.Get('h_fj1MSD_tt_matched')
hunmatched = ftemplate.Get('h_fj1MSD_tt_unmatched')
hW = ftemplate.Get('h_fj1MSD_W')

mass = root.RooRealVar("m","m_{SD} [GeV]",40,450)

normbound=0.2

norm_ = hunmatched.Integral()
ttu_norm = root.RooRealVar('ttu_norm','ttu_norm',norm_,(1-normbound)*norm_,(1+normbound)*norm_)
dh_ttu = root.RooDataHist('dh_ttu','dh_ttu',root.RooArgList(mass),hunmatched)
ttu = root.RooHistPdf('ttu','ttu',root.RooArgSet(mass),dh_ttu)

norm_ = hmatched.Integral()
ttm_norm = root.RooRealVar('ttm_norm','ttm_norm',norm_,(1-normbound)*norm_,(1+normbound)*norm_)
dh_ttm = root.RooDataHist('dh_ttm','dh_ttm',root.RooArgList(mass),hmatched)
ttm = root.RooHistPdf('ttm','ttm',root.RooArgSet(mass),dh_ttm)

norm_ = hW.Integral()
W_norm = root.RooRealVar('W_norm','W_norm',norm_,(1-normbound)*norm_,(1+normbound)*norm_)
dh_W = root.RooDataHist('dh_W','dh_W',root.RooArgList(mass),hW)
W = root.RooHistPdf('W','W',root.RooArgSet(mass),dh_W)

model = root.RooAddPdf('model','model',root.RooArgList(ttm,ttu,W),root.RooArgList(ttm_norm,ttu_norm,W_norm))

dh_data = root.RooDataHist('dh_data','dh_data',root.RooArgList(mass),hdata)
# pdf_data = root.RooHistPdf('datapdf','datapdf',root.RooArgSet(mass),dh_data,0)

model.fitTo(dh_data)

frame = mass.frame()
dh_data.plotOn(frame)
hdata.Scale((ttu_norm.getVal()+ttm_norm.getVal()+W_norm.getVal())/hdata.Integral())
model.plotOn(frame)
#for c,l in [('ttm',root.RooFit.LineColor(root.kRed))]:
for c,l in [('ttm',root.RooFit.LineColor(root.kOrange)),('ttu',root.RooFit.LineColor(6)),('W',root.RooFit.LineColor(8))]:
  model.plotOn(frame,root.RooFit.Components(c),l)

c = root.TCanvas()
frame.Draw()
for h in [hmatched,hunmatched,hW]:
  h.SetLineWidth(1)
  h.SetLineStyle(2)
hmatched.SetLineColor(root.kOrange); hmatched.Draw('samehist')
hunmatched.SetLineColor(6); hunmatched.Draw('samehist')
hW.SetLineColor(8); hW.Draw('samehist')

for ext in ['pdf','png']:
  c.SaveAs('~/public_html/figs/toptagging/datavalidation/v7/templates/simplefitMC.'+ext)
