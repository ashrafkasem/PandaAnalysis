#!/usr/bin/env python

from sys import argv
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
Load('Drawers','HistogramDrawer')

plot = {}
for iC in [0,1]:
  plot[iC] = root.HistogramDrawer()
  plot[iC].Ratio(1)
  plot[iC].FixRatio(.5)
  plot[iC].DrawMCErrors(False)
  plot[iC].Stack(True)
  plot[iC].DrawEmpty(True)
  plot[iC].SetTDRStyle()
  plot[iC].AddCMSLabel()
  plot[iC].SetLumi(12.9); plot[iC].AddLumiLabel(True)
  plot[iC].InitLegend()

# masscorr = 'L2L3'
masscorr = ''

hprong = {}; dhprong = {}; pdfprong = {}; norm = {}; smeared = {}; smear = {}; mu = {}; hdata = {}; dh_data={}
mass = root.RooRealVar("m","m_{SD} [GeV]",50,450)

ftemplate = {
      'pass' : root.TFile('~/public_html/figs/toptagging/datavalidation/v8/templates/tag_top_ecfv8_bdt_pass_hists.root'),
      'fail' : root.TFile('~/public_html/figs/toptagging/datavalidation/v8/templates/tag_top_ecfv8_bdt_fail_hists.root')
    }

# get histograms
hdata[1] = ftemplate['pass'].Get('h_fj1MSD%s_Data'%masscorr)
dh_data[1] = root.RooDataHist('dh_data1','dh_data1',root.RooArgList(mass),hdata[1])
hprong[(3,1)] = ftemplate['pass'].Get('h_fj1MSD%s_3-prong'%masscorr)
hprong[(2,1)] = ftemplate['pass'].Get('h_fj1MSD%s_2-prong'%masscorr)
hprong[(1,1)] = ftemplate['pass'].Get('h_fj1MSD%s_1-prong'%masscorr)

hdata[0] = ftemplate['fail'].Get('h_fj1MSD%s_Data'%masscorr)
dh_data[0] = root.RooDataHist('dh_data0','dh_data0',root.RooArgList(mass),hdata[0])
hprong[(3,0)] = ftemplate['fail'].Get('h_fj1MSD%s_3-prong'%masscorr)
hprong[(2,0)] = ftemplate['fail'].Get('h_fj1MSD%s_2-prong'%masscorr)
hprong[(1,0)] = ftemplate['fail'].Get('h_fj1MSD%s_1-prong'%masscorr)

# build pdfs
for iC in [0,1]:
  for iP in xrange(1,4):
    cat = (iP,iC)
    dhprong[cat] = root.RooDataHist('dh%i%i'%cat,'dh%i%i'%cat,root.RooArgList(mass),hprong[cat]) 
    pdfprong[cat] = root.RooHistPdf('pdf%i%i'%cat,'pdf%i%i'%cat,root.RooArgSet(mass),dhprong[cat]) 
    norm_ = hprong[cat].Integral()
    norm[cat] = root.RooRealVar('norm%i%i'%cat,'norm%i%i'%cat,norm_,0.01*norm_,100*norm_)

# smear pdfs
dummy = root.RooRealVar('dummy','dummy',0.1,0.1,10)
for iP in [1,2,3]:
  mu[iP] = root.RooRealVar('mu%i'%iP,'mu%i'%iP,10,-25,25)
  smear[iP] = root.RooGaussian('gauss%i'%iP,'gauss%i'%iP,mass,mu[iP],dummy)
  for iC in [0,1]:
    cat = (iP,iC)
    if iP<3:
      smeared[cat] = pdfprong[cat]
    else:
      smeared[cat] = root.RooFFTConvPdf('conv%i%i'%cat,'conv%i%i'%cat,mass,pdfprong[cat],smear[iP])

model = {}
nsigtotal = root.RooFormulaVar('nsigtotal','norm30+norm31',root.RooArgList(norm[(3,0)],norm[(3,1)]))
eff_ = norm[(3,1)].getVal()/(norm[(3,1)].getVal()+norm[(3,0)].getVal())
eff = root.RooRealVar('eff','eff',eff_,0.5*eff_,2*eff_)
normsig = {
      0 : root.RooFormulaVar('nsigfail','(1.0-eff)*nsigtotal',root.RooArgList(eff,nsigtotal)),
      1 : root.RooFormulaVar('nsigpass','eff*nsigtotal',root.RooArgList(eff,nsigtotal)),
    }
for iC in [0,1]:
  model[iC] = root.RooAddPdf('model%i'%iC,'model%i'%iC,
                            root.RooArgList(*[smeared[(x,iC)] for x in [1,2,3]]),
                            root.RooArgList(norm[(1,iC)],norm[(2,iC)],normsig[iC]))

# build simultaneous fit
sample = root.RooCategory('sample','')
sample.defineType('pass',1)
sample.defineType('fail',2)

datacomb = root.RooDataHist('datacomb','datacomb',root.RooArgList(mass),
                            root.RooFit.Index(sample),
                            root.RooFit.Import('pass',dh_data[1]),
                            root.RooFit.Import('fail',dh_data[0]))

simult = root.RooSimultaneous('simult','simult',sample)
simult.addPdf(model[1],'pass')
simult.addPdf(model[0],'fail')

# fit!
fitresult = simult.fitTo(datacomb,
                          root.RooFit.Extended(),
                          root.RooFit.Strategy(2),
                          root.RooFit.Minos(root.RooArgSet(eff)),
                          root.RooFit.NumCPU(4),
                          root.RooFit.Save())


# dump the efficiencies
pcat=(3,1); fcat=(3,0)
masslo=110; masshi =210 # to get on binedges for now
effMass_ = hprong[pcat].Integral(hprong[pcat].FindBin(masslo),hprong[pcat].FindBin(masshi)-1)/hprong[pcat].Integral() # prefit efficiency of the mass cut on the pass distribution
mass.setRange('MASSWINDOW',masslo,masshi)
massint = smeared[(3,1)].createIntegral(root.RooArgSet(mass),root.RooArgSet(mass),'MASSWINDOW')
effMass = massint.getVal(); effMassErr = massint.getPropagatedError(fitresult)

print 'Tagging cut:'
print '\tPre-fit efficiency was %f'%(eff_)
print '\tPost-fit efficiency is %f +%.5g -%.5g'%(eff.getVal(),abs(eff.getErrorHi()),abs(eff.getErrorLo()))
print 'Tagging+mass cut:'
print '\tPre-fit mass efficiency was %f'%(effMass_)
print '\tPost-fit mass efficiency is %f +/-%.5g'%(effMass,effMassErr)
print '\tPre-fit mass+tag efficiency was %f'%(effMass_*eff_)
print '\tPost-fit mass+tag efficiency is %f +%.5g -%.5g'%(effMass*eff.getVal(),effMass*abs(eff.getErrorHi()),effMass*abs(eff.getErrorLo()))

# make nice plots

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

for iC in [0,1]:
  for iP in [1,2,3]:
    cat = (iP,iC)
    h = smeared[cat].createHistogram('h%i%i'%cat,mass,root.RooFit.Binning(40))
    h.SetLineWidth(3)
    h.SetLineStyle(1)
    h.SetLineColor(colors[iP])
    h.Scale(norm[cat].getVal()/h.Integral())
    plot[iC].AddAdditional(h,'hist',labels[iP])
    
  hprefit = hprong[(1,0)].Clone('prefit')
  hprefit.Reset()
  for jP in [1,2,3]:
    hprefit.Add(hprong[(jP,iC)])
  hprefit.SetLineWidth(2)
  hprefit.SetLineStyle(2)
  hprefit.SetLineColor(root.kBlue+2)
  plot[iC].AddAdditional(hprefit,'hist','Pre-fit')

  for iP in [1,2,3]:
    cat = (iP,iC)
    hprong[cat].SetLineColor(colors[iP])
    hprong[cat].SetLineWidth(2)
    hprong[cat].SetLineStyle(2)
    plot[iC].AddAdditional(hprong[cat],'hist')

  hdata[iC].SetLineColor(root.kBlack)
  hdata[iC].SetMarkerStyle(20);
  plot[iC].AddHistogram(hdata[iC],'Data',root.kData)

  hmodel_ = model[iC].createHistogram('hmodel%i'%iC,mass,root.RooFit.Binning(40))
  hmodel = root.TH1D(); hmodel_.Copy(hmodel)
  hmodel.SetLineWidth(3);
  hmodel.SetLineColor(root.kBlue+10)
  if masscorr=='':
    hmodel.GetXaxis().SetTitle('fatjet m_{SD} [GeV]')
  else:
    hmodel.GetXaxis().SetTitle('L2L3-corr fatjet m_{SD} [GeV]')
  hmodel.GetYaxis().SetTitle('Events/10 GeV')
  hmodel.Scale(sum([norm[(x,iC)].getVal() for x in [1,2,3]])/hmodel.Integral())
  plot[iC].AddHistogram(hmodel,'Post-fit',root.kExtra5)
  plot[iC].AddAdditional(hmodel,'hist')

plot[1].AddPlotLabel('Pass category',.18,.77,False,42,.04)
plot[0].AddPlotLabel('Fail category',.18,.77,False,42,.04)
plot[1].Draw('~/public_html/figs/toptagging/datavalidation/v8/templates/','pass%s'%masscorr)
plot[0].Draw('~/public_html/figs/toptagging/datavalidation/v8/templates/','fail%s'%masscorr)

# # mass.setVal(175)
# # print pdfprong[iP].getVal(root.RooArgSet(mass))
# # print smear[iP].getVal(root.RooArgSet(mass))
# # print smeared[iP].getVal(root.RooArgSet(mass))
