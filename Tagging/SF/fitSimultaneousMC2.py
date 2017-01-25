#!/usr/bin/env python

import argparse
from sys import argv,exit

parser = argparse.ArgumentParser(description='fit stuff')
parser.add_argument('--indir',metavar='indir',type=str)
parser.add_argument('--jesr',metavar='jesr',type=str,default='central')
args = parser.parse_args()
basedir = args.indir+'/'+args.jesr+'/'
argv=[]

from math import sqrt
import ROOT as root
from PandaCore.Tools.Load import *
Load('Drawers','HistogramDrawer')

def imp(w_):
	return getattr(w_,'import')

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
	plot[iC].SetLumi(36.6)
	plot[iC].AddLumiLabel(True)
	plot[iC].InitLegend()

# masscorr = 'L2L3'
masscorr = ''
postfix = {
		'central' : '',
		'scaleUp' : 'ScaleUp',
		'scaleDown' : 'ScaleDown',
		'smeared' : 'Smeared'
		}
masscorr += postfix[args.jesr]

hprong = {}; dhprong = {}; pdfprong = {}; norm = {}; smeared = {}; smear = {}; mu = {}; hdata = {}; dh_data={}
mcnorms = {}; mcerrs = {}
mass = root.RooRealVar("m","m_{SD} [GeV]",50,450)

ftemplate = {
			'pass' : root.TFile(basedir+'top_ecf_bdt_pass_hists.root'),
			'fail' : root.TFile(basedir+'top_ecf_bdt_fail_hists.root')
		}

# get histograms
hdata[1] = ftemplate['pass'].Get('h_fj1MSD%s_Data'%masscorr)
hprong[(3,1)] = ftemplate['pass'].Get('h_fj1MSD%s_3-prong'%masscorr)
hprong[(2,1)] = ftemplate['pass'].Get('h_fj1MSD%s_2-prong'%masscorr)
hprong[(1,1)] = ftemplate['pass'].Get('h_fj1MSD%s_1-prong'%masscorr)
hprong[(0,1)] = ftemplate['pass'].Get('h_fj1MSD%s_bkg'%masscorr)
#hprong[(0,1)] = hprong[(1,1)].Clone(); 
#hprong[(0,1)].Add(hprong[(2,1)])

hdata[0] = ftemplate['fail'].Get('h_fj1MSD%s_Data'%masscorr)
hprong[(3,0)] = ftemplate['fail'].Get('h_fj1MSD%s_3-prong'%masscorr)
hprong[(2,0)] = ftemplate['fail'].Get('h_fj1MSD%s_2-prong'%masscorr)
hprong[(1,0)] = ftemplate['fail'].Get('h_fj1MSD%s_1-prong'%masscorr)
hprong[(0,0)] = ftemplate['fail'].Get('h_fj1MSD%s_bkg'%masscorr)
#hprong[(0,0)] = hprong[(1,0)].Clone(); hprong[(0,0)].Add(hprong[(2,0)])

NBINS = hdata[1].GetNbinsX()

for iB in xrange(1,NBINS+1):
	# hack to add the MC unc as a data unc
	for iC in [0,1]:
		err_0 = hprong[(0,iC)].GetBinError(iB)
		err_3 = hprong[(3,iC)].GetBinError(iB)
		err_d = hdata[iC].GetBinError(iB)
		err = sqrt(pow(err_0,2)+pow(err_3,2)+pow(err_d,2))
		hdata[iC].SetBinError(iB,err)

dh_data[1] = root.RooDataHist('dh_data1','dh_data1',root.RooArgList(mass),hdata[1])
dh_data[0] = root.RooDataHist('dh_data0','dh_data0',root.RooArgList(mass),hdata[0])

# build pdfs
for iC in [0,1]:
	for iP in [3,0]:
		cat = (iP,iC)
		print '\033[0;41m                                                              \033[0m'
		print cat, hprong[cat].Integral()
		dhprong[cat] = root.RooDataHist('dh%i%i'%cat,'dh%i%i'%cat,root.RooArgList(mass),hprong[cat]) 
		pdfprong[cat] = root.RooHistPdf('pdf%i%i'%cat,'pdf%i%i'%cat,root.RooArgSet(mass),dhprong[cat]) 
		norm_ = hprong[cat].Integral()
		#norm[cat] = root.RooRealVar('norm%i%i'%cat,'norm%i%i'%cat,norm_,0.01*norm_,100*norm_)
		norm[cat] = root.RooRealVar('norm%i%i'%cat,'norm%i%i'%cat,norm_,0.7*norm_,1.3*norm_)
		mcnorms[cat] = norm_
		err_ = 0
		for iB in xrange(1,hprong[cat].GetNbinsX()+1):
			err_ += pow(hprong[cat].GetBinError(iB),2)
		mcerrs[cat] = sqrt(err_)

for cat in sorted(mcnorms):
	print cat,mcnorms[cat]
print '\033[0;41m                                                              \033[0m'

mass.setBins(2000,'cache')
# smear pdfs with a single gaussian for pass/fail
jesr_sigma = {x:root.RooRealVar('jesr_sigma%i'%x,'jesr_sigma%i'%x,0.5,0.05,10) for x in [0,1,3]}
xmus_all = { 
		'central'   : {0:(-2,-10,0),1:(0,-10,10),3:(0,-10,10)},
		'smeared'   : {0:(-2,-10,0),1:(0,-10,10),3:(0,-10,10)},
		'scaleUp'   : {0:(-2,-10,0),1:(0,-10,10),3:(0,-10,10)},
		'scaleDown' : {0:(-2,-10,0),1:(0,-10,10),3:(0,-10,10)},
		}
xmus = xmus_all[args.jesr]
jesr_mu = {x:root.RooRealVar('jesr_mu%i'%x,'jesr_mu%i'%x,xmus[x][0],xmus[x][1],xmus[x][2]) for x in [0,1,3]}
jesr = {x:root.RooGaussian('jesr%i'%x,'jesr%i'%x,mass,jesr_mu[x],jesr_sigma[x]) for x in [0,1,3]}
for iP in [3,0]:
	for iC in [0,1]:
		cat = (iP,iC)
		if iP==3 and iC==0:
# 	 if False:
			smeared[cat] = root.RooFFTConvPdf('conv%i%i'%cat,'conv%i%i'%cat,mass,pdfprong[cat],jesr[iC])
		else:
			smeared[cat] = pdfprong[cat]

model = {}
nsigtotal = root.RooFormulaVar('nsigtotal','norm30+norm31',root.RooArgList(norm[(3,0)],norm[(3,1)]))
# eff = root.RooFormulaVar('eff','norm31/(norm30+norm31)',root.RooArgList(norm[(3,1)],norm[(3,0)]))
# eff_ = norm[(3,1)].getVal()/(norm[(3,1)].getVal()+norm[(3,0)].getVal())

def calcEffAndErr(p,perr,f,ferr):
	eff_ = p/(p+f)
	err_ = pow( perr * f / pow(p+f,2) , 2 )
	err_ += pow( ferr * p / pow(p+f,2) , 2 )
	err_ = sqrt(err_)
	return eff_,err_

eff_,err_ = calcEffAndErr(mcnorms[(3,1)],mcerrs[(3,1)],mcnorms[(3,0)],mcnorms[(3,1)])
eff_ = mcnorms[(3,1)]/(mcnorms[(3,1)]+mcnorms[(3,0)])
eff = root.RooRealVar('eff','eff',eff_,0.66*eff_,1.5*eff_)
normsig = {
	    0 : root.RooFormulaVar('nsigfail','(1.0-eff)*nsigtotal',root.RooArgList(eff,nsigtotal)),
	    1 : root.RooFormulaVar('nsigpass','eff*nsigtotal',root.RooArgList(eff,nsigtotal)),
	  }
for iC in [0,1]:
	model[iC] = root.RooAddPdf('model%i'%iC,'model%i'%iC,
#                            root.RooArgList(*[pdfprong[(x,iC)] for x in [0,3]]),
	                          root.RooArgList(*[smeared[(x,iC)] for x in [0,3]]),
#                            root.RooArgList(norm[(0,iC)],norm[(3,iC)]))
	                          root.RooArgList(norm[(0,iC)],normsig[iC]))

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
	                        #root.RooFit.Minos(root.RooArgSet(norm[(3,0)],norm[(3,1)])),
	                        root.RooFit.Minos(root.RooArgSet(eff)),
	                        root.RooFit.NumCPU(4),
	                        root.RooFit.Save())


# dump the efficiencies
pcat=(3,1); fcat=(3,0)
masslo=110; masshi =210 
effMass_ = hprong[pcat].Integral(hprong[pcat].FindBin(masslo),hprong[pcat].FindBin(masshi)-1)/hprong[pcat].Integral() # prefit efficiency of the mass cut on the pass distribution
mass.setRange('MASSWINDOW',masslo,masshi)
massint = smeared[(3,1)].createIntegral(root.RooArgSet(mass),root.RooArgSet(mass),'MASSWINDOW')
effMass = massint.getVal(); effMassErr = massint.getPropagatedError(fitresult)

# make nice plots

colors = {
	1:8,
	2:6,
	3:root.kOrange-3,
	0:root.kCyan+3,
}
labels = {
	1:'W+jets',
	2:'Unmatched top/VV',
	3:'Matched top',
	0:'Background',
}

for iC in [0,1]:
	for iP in [0,3]:
	  cat = (iP,iC)
	  h = smeared[cat].createHistogram('h%i%i'%cat,mass,root.RooFit.Binning(NBINS))
	  h.SetLineWidth(3)
	  h.SetLineStyle(1)
	  h.SetLineColor(colors[iP])
	  h.Scale(norm[cat].getVal()/h.Integral())
	  plot[iC].AddAdditional(h,'hist',labels[iP])
	  
	hprefit = hprong[(1,0)].Clone('prefit')
	for iB in xrange(1,hprefit.GetNbinsX()+1):
	  hprefit.SetBinContent(iB,0)
	hprefit.Reset()
	for jP in [3,0]:
	  hprefit.Add(hprong[(jP,iC)])
	hprefit.SetLineWidth(2)
	hprefit.SetLineStyle(2)
	hprefit.SetLineColor(root.kBlue+2)
	hprefit.SetFillStyle(0)
	plot[iC].AddAdditional(hprefit,'hist','Pre-fit')

	for iP in [0,3]:
	  cat = (iP,iC)
	  hprong[cat].SetLineColor(colors[iP])
	  hprong[cat].SetFillStyle(0)
	  hprong[cat].SetLineWidth(2)
	  hprong[cat].SetLineStyle(2)
	  plot[iC].AddAdditional(hprong[cat],'hist')

	hdata[iC].SetLineColor(root.kBlack)
	hdata[iC].SetMarkerStyle(20);
	plot[iC].AddHistogram(hdata[iC],'Data',root.kData)

	hmodel_ = model[iC].createHistogram('hmodel%i'%iC,mass,root.RooFit.Binning(NBINS))
	hmodel = root.TH1D(); hmodel_.Copy(hmodel)
	hmodel.SetLineWidth(3);
	hmodel.SetLineColor(root.kBlue+10)
	if masscorr=='':
	  hmodel.GetXaxis().SetTitle('fatjet m_{SD} [GeV]')
	else:
	  hmodel.GetXaxis().SetTitle('L2L3-corr fatjet m_{SD} [GeV]')
	hmodel.GetYaxis().SetTitle('Events/10 GeV')
	hmodel.Scale(sum([norm[(x,iC)].getVal() for x in [0,3]])/hmodel.Integral())
	hmodel.SetFillStyle(0)
	plot[iC].AddHistogram(hmodel,'Post-fit',root.kExtra5)
	plot[iC].AddAdditional(hmodel,'hist')

	eff_post = eff.getVal()
	#eff_posthi = eff.getPropagatedError(fitresult)
	#eff_postlo = eff.getPropagatedError(fitresult)
	eff_posthi = eff.getErrorHi()
	eff_postlo = abs(eff.getErrorLo())
	eff_err = max(eff_posthi,eff_postlo)

	plot[iC].AddPlotLabel('#varepsilon_{tag}^{Data} = %.3g#pm%.3g'%(eff_post,eff_err),
	                      .6,.47,False,42,.04)
	#plot[iC].AddPlotLabel('#varepsilon_{tag}^{MC} = %.3g^{+%.2g}_{-%.2g}'%(eff_,err_,err_),
	plot[iC].AddPlotLabel('#varepsilon_{tag}^{MC} = %.3g'%(eff_),
	                      .6,.37,False,42,.04)
	plot[iC].AddPlotLabel('#varepsilon_{tag+mSD}^{Data} = %.3g#pm%.3g'%(effMass*eff_post,effMass*eff_err),
	                      .6,.27,False,42,.04)
	#plot[iC].AddPlotLabel('#varepsilon_{tag+mSD}^{MC} = %.3g^{+%.2g}_{-%.2g}'%(eff_*effMass_,err_*effMass_,err_*effMass_),
	plot[iC].AddPlotLabel('#varepsilon_{tag+mSD}^{MC} = %.3g'%(eff_*effMass_),
	                      .6,.17,False,42,.04)

plot[1].AddPlotLabel('Pass category',.18,.77,False,42,.05)
plot[0].AddPlotLabel('Fail category',.18,.77,False,42,.05)
plot[1].Draw(basedir,'pass%s'%masscorr)
plot[0].Draw(basedir,'fail%s'%masscorr)

# save output
#rcn = root.RooFit.Silence
w = root.RooWorkspace('w','workspace')
w.imp = imp(w)
w.imp(mass)
for x in [nsigtotal,eff,sample,datacomb,simult]:
	w.imp(x,root.RooFit.Silence())
for iC in [0,1]:
	w.imp(normsig[iC],root.RooFit.Silence())
	w.imp(dh_data[iC],root.RooFit.Silence())
	w.imp(model[iC],root.RooFit.Silence())
	for iP in [0,3]:
	  cat = (iP,iC)
#    w.imp(smear[iP]); w.imp(mu[iP])
	  w.imp(dhprong[cat],root.RooFit.Silence())
	  w.imp(pdfprong[cat],root.RooFit.Silence())
	  w.imp(norm[cat],root.RooFit.Silence())
#    w.imp(smeared[cat])
w.writeToFile(basedir+'wspace.root')

print 'Tagging cut:'
print '\tPre-fit efficiency was %f'%(eff_)
print '\tPost-fit efficiency is %f +%.5g -%.5g'%(eff_post,eff_posthi,eff_postlo)
print 'Tagging+mass cut:'
print '\tPre-fit mass efficiency was %f'%(effMass_)
print '\tPost-fit mass efficiency is %f +/-%.5g'%(effMass,effMassErr)
print '\tPre-fit mass+tag efficiency was %f +/%.2g'%(effMass_*eff_,effMass_*err_)
print '\tPost-fit mass+tag efficiency is %f +%.5g -%.5g'%(effMass*eff_post,effMass*eff_posthi,effMass*eff_postlo)
