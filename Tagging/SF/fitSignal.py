#!/usr/bin/env python

from sys import argv
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
Load('Drawers','HistogramDrawer')

plot = root.HistogramDrawer()
plot.Ratio(1)
plot.FixRatio(.5)
plot.DrawMCErrors(False)
plot.Stack(True)
plot.DrawEmpty(True)
plot.SetTDRStyle()
plot.AddCMSLabel()
plot.SetLumi(12.9); plot.AddLumiLabel(True)
plot.InitLegend()

#masscorr = 'L2L3'
masscorr = ''

ftemplate = root.TFile('~/public_html/figs/toptagging/datavalidation/v8/templates/tag_top_ecfv8_bdt_pass_hists.root')
hdata = ftemplate.Get('h_fj1MSD%s_Data'%masscorr)
hprong = {}; dhprong = {}; pdfprong = {}; norm = {}; smeared = {}; smear = {}; mu = {}
hprong[3] = ftemplate.Get('h_fj1MSD%s_3-prong'%masscorr)
hprong[2] = ftemplate.Get('h_fj1MSD%s_2-prong'%masscorr)
hprong[1] = ftemplate.Get('h_fj1MSD%s_1-prong'%masscorr)

mass = root.RooRealVar("m","m_{SD} [GeV]",50,450)

for iP in xrange(1,4):
  dhprong[iP] = root.RooDataHist('dh%i'%iP,'dh%i'%iP,root.RooArgList(mass),hprong[iP]) 
  pdfprong[iP] = root.RooHistPdf('pdf%i'%iP,'pdf%i'%iP,root.RooArgSet(mass),dhprong[iP]) 
  norm_ = hprong[iP].Integral()
  norm[iP] = root.RooRealVar('normm%i'%iP,'norm%i'%iP,norm_,0.1*norm_,10*norm_)

dh_data = root.RooDataHist('dh_data','dh_data',root.RooArgList(mass),hdata)

dummy = root.RooRealVar('dummy','dummy',0.1,0.1,10)
for iP in [1,2,3]:
  if iP<3: 
    smeared[iP] = pdfprong[iP]
  else:
    mu[iP] = root.RooRealVar('mu%i'%iP,'mu%i'%iP,10,-25,25)
    smear[iP] = root.RooGaussian('gauss%i'%iP,'gauss%i'%iP,mass,mu[iP],dummy)
    smeared[iP] = root.RooFFTConvPdf('conv%i'%iP,'conv%i'%iP,mass,pdfprong[iP],smear[iP])

model = root.RooAddPdf('model','model',root.RooArgList(*[smeared[x] for x in [1,2,3]]),root.RooArgList(*[norm[x] for x in [1,2,3]]))

model.fitTo(dh_data)


c = root.TCanvas()

colors = {
  1:8,
  2:6,
  3:root.kOrange
}
labels = {
  1:'W+jets',
  2:'Unmatched top/VV',
  3:'Matched top',
}
for iP in [1,2,3]:
  h = smeared[iP].createHistogram('h%i'%iP,mass,root.RooFit.Binning(40))
  h.SetLineWidth(3)
  h.SetLineStyle(1)
  h.SetLineColor(colors[iP])
  h.Scale(norm[iP].getVal()/h.Integral())
  plot.AddAdditional(h,'hist',labels[iP])

hprong[3].SetLineColor(colors[3]); 
hprong[2].SetLineColor(colors[2]); 
hprong[1].SetLineColor(colors[1]); 
for h in [hprong[x] for x in [3,2,1]]:
  h.SetLineWidth(2)
  h.SetLineStyle(2)
  plot.AddAdditional(h,'hist')

hdata.SetLineColor(root.kBlack)
hdata.SetMarkerStyle(20);
plot.AddHistogram(hdata,'Data',root.kData)

hmodel_ = model.createHistogram('hmodel',mass,root.RooFit.Binning(40))
hmodel = root.TH1D(); hmodel_.Copy(hmodel)
hmodel.SetLineWidth(3);
hmodel.SetLineColor(root.kBlue+10)
if masscorr=='':
  hmodel.GetXaxis().SetTitle('fatjet m_{SD} [GeV]')
else:
  hmodel.GetXaxis().SetTitle('L2L3-corr fatjet m_{SD} [GeV]')
hmodel.GetYaxis().SetTitle('Events/10.25 GeV')
hmodel.Scale(sum([x.getVal() for x in norm.values()])/hmodel.Integral())
plot.AddHistogram(hmodel,'Post-fit',root.kExtra5)
plot.AddAdditional(hmodel,'hist')

hprefit = hprong[1].Clone('prefit')
hprefit.Reset()
for iP in [1,2,3]:
  hprefit.Add(hprong[iP])
hprefit.SetLineWidth(2)
hprefit.SetLineStyle(2)
hprefit.SetLineColor(root.kBlue+2)
plot.AddAdditional(hprefit,'hist','Pre-fit')

plot.Draw('~/public_html/figs/toptagging/datavalidation/v8/templates/','simplefit%s'%masscorr)

# mass.setVal(175)
# print pdfprong[iP].getVal(root.RooArgSet(mass))
# print smear[iP].getVal(root.RooArgSet(mass))
# print smeared[iP].getVal(root.RooArgSet(mass))
