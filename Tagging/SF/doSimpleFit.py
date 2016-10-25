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

normbound=0.5

norm_ = hunmatched.Integral()
ttu_norm = root.RooRealVar('ttu_norm','ttu_norm',norm_,(1-normbound)*norm_,(1+normbound)*norm_)
ttu_mu = root.RooRealVar('ttu_mu','mean ttu',92,92,92)
ttu_sigma = root.RooRealVar('ttu_sigma','sigma ttu',25,25,25)
ttu = root.RooLandau('ttu','ttu',mass,ttu_mu,ttu_sigma)

norm_ = hmatched.Integral()
ttm_norm = root.RooRealVar('ttm_norm','ttm_norm',norm_,(1-normbound)*norm_,(1+normbound)*norm_)
ttm_mu = root.RooRealVar('ttm_mu','mean ttm',162,140,200)
ttm_sigma = root.RooRealVar('ttm_sigma','sigma ttm',65,63,67)
ttm = root.RooBreitWigner('ttm','ttm',mass,ttm_mu,ttm_sigma)

norm_ = hW.Integral()
W_norm = root.RooRealVar('W_norm','W_norm',norm_,(1-normbound)*norm_,(1+normbound)*norm_)
W_a = root.RooRealVar('W_a','W_a',-0.03,-1,-0.00001)
W_mu = root.RooRealVar('W_mu','W_mu',0,0,1000)
W_sigma = root.RooRealVar('W_sigma','W_sigma',50,0,1000)
#W = root.RooLognormal('W','W',mass,W_mu,W_sigma)
W = root.RooExponential('W','W',mass,W_a)

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
hmatched.SetLineWidth(3); hmatched.SetLineColor(root.kOrange); hmatched.Draw('samehist')
hunmatched.SetLineWidth(3); hunmatched.SetLineColor(6); hunmatched.Draw('samehist')
hW.SetLineWidth(3); hW.SetLineColor(8); hW.Draw('samehist')

for ext in ['pdf','png']:
  c.SaveAs('~/public_html/figs/toptagging/datavalidation/v7/templates/simplefit.'+ext)
