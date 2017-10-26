#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse

basedir = getenv('PANDA_FLATDIR')

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--indir',metavar='indir',type=str,default=basedir)
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--wrt',metavar='wrt',type=str,default=None)
parser.add_argument('--sel',metavar='sel',type=str,default='inc') # inc, tag, mistag
args = parser.parse_args()

figsdir = args.outdir
basedir = args.indir
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
Load('HistogramDrawer')

lumi=35800
plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.SetLumi(lumi/1000.)
plot.AddCMSLabel()
plot.AddLumiLabel()
plot.SetAutoRange(False)
plot.InitLegend()
plot.AddPlotLabel('110 < m_{SD} < 210 GeV',.18,.77,False,42,.04)
# root.gStyle.SetPadRightMargin(0.15)
root.gStyle.SetOptStat(0)

nBins=20

tmc = root.TChain('events')
tmc.AddFile(basedir+'/TTbar.root')
tmc.AddFile(basedir+'/WJets.root')
fdata = root.TFile(basedir+'/MET.root')
tdata = fdata.Get('events')

cut = 'nFatjet==1 && fj1Pt>250 && nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && pfUWmag>250 && fj1MSD>110 && fj1MSD<210'
weight = '%f*normalizedWeight*sf_pu*sf_lepID*sf_lepIso*sf_ewkV*sf_qcdV*sf_tt'%lumi
label = ''
if args.sel=='tag':
  cut = 'nFatjet==1 && fj1Pt>250 && nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && pfUWmag>250 && fj1MSD>110 && fj1MSD<210 && fj1MaxCSV>0.54 && isojetNBtags==1'
  weight = '%f*normalizedWeight*sf_pu*sf_lepID*sf_lepIso*sf_ewkV*sf_qcdV*sf_btag1*sf_sjbtag1*sf_tt'%lumi
  label = 'tag_'
elif args.sel=='mistag':
  cut = 'nFatjet==1 && fj1Pt>250 && nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && pfUWmag>250 && fj1MSD>110 && fj1MSD<210 && fj1MaxCSV<0.54 && isojetNBtags==0'
  weight = '%f*normalizedWeight*sf_pu*sf_lepID*sf_lepIso*sf_ewkV*sf_qcdV*sf_btag0*sf_sjbtag0*sf_tt'%lumi
  label = 'mistag_'

if args.wrt=='npv':
  wrtvar='npv'
  wrtlo=5
  wrthi=35
  wrtbins=30
  xaxis='NPV'
elif args.wrt=='pt':
  wrtvar='fj1Pt'
  wrtlo=250
  wrthi=750
  wrtbins=40
  xaxis='Fatjet p_{T} [GeV]'
label += args.wrt+'_'

taggervar='top_ecf_bdt'
taggerlo=-1.2
taggerhi=1.
taggerbins=100
yaxis = '#LTBDT#GT'
hmc = root.TH2D('hmc','hmc',wrtbins,wrtlo,wrthi,taggerbins,taggerlo,taggerhi)
hdata = root.TH2D('hdata','hdata',wrtbins,wrtlo,wrthi,taggerbins,taggerlo,taggerhi)
tmc.Draw('%s:%s>>hmc'%(taggervar,wrtvar),tTIMES(weight,cut),'colz')
tdata.Draw('%s:%s>>hdata'%(taggervar,wrtvar),tAND('(trigger&1)!=0',cut),'colz')

pmc = hmc.ProfileX()
pdata = hdata.ProfileX()
pdata.SetMinimum(-0.6); pdata.SetMaximum(.6)
pdata.GetXaxis().SetTitle(xaxis)
pdata.GetYaxis().SetTitle(yaxis)

pmc.SetFillStyle(3001)
pmc.SetFillColorAlpha(root.kRed,1)
pmc.SetLineWidth(3); pmc.SetLineColor(root.kRed)
plot.AddAdditional(pmc,'e2','MC')

pdata.SetFillStyle(0)
pdata.SetLineWidth(2); pdata.SetLineColor(root.kBlack)
plot.AddHistogram(pdata,'Data',root.kData)
plot.AddAdditional(pdata,'elp')

plot.Draw(args.outdir+'/',label+taggervar)
