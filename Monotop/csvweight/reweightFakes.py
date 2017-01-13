#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
parser = argparse.ArgumentParser(description='reweight CSV shape for fakes')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
args = parser.parse_args()
lumi = 36.56
blind=True
linear=False
sname = argv[0]

argv=[]
import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import *
import PandaCore.Tools.Functions
Load('Drawers','HistogramDrawer')

fin = root.TFile(args.outdir+'/photon_hists.root')

hdata = fin.Get('h_isojet1CSV_Data')
hglf = fin.Get('h_isojet1CSV_#gamma+jets (LF)')
hghf = fin.Get('h_isojet1CSV_#gamma+jets (HF)')
hqcd = fin.Get('h_isojet1CSV_QCD')
hfake = hglf.Clone('h_isojet1CSV_bkg')
hfake.Add(hqcd)
hfake.SetFillStyle(0)

hdata.Add(hghf,-1) # this is what we "estimate" to be fakes
hdata.Scale(1./hdata.Integral())
hfake.Scale(1./hfake.Integral())

hratio = hdata.Clone('hratio')
hratio.Divide(hfake)
hratio.GetYaxis().SetTitle('Data/MC')
hratio.SetMaximum(2.5)
hratio.SetMinimum(0)
hratioerr = hratio.Clone()
hratioerr.SetFillStyle(3005)
hratioerr.SetFillColorAlpha(root.kBlack,0.5)
hratioerr.SetLineWidth(0)
hratioerr.SetMarkerStyle(0)

### LOAD PLOTTING UTILITY ###
plot = root.HistogramDrawer()
plot.DrawEmpty(True)
plot.SetLumi(lumi)
plot.SetTDRStyle()
plot.AddCMSLabel()
plot.AddLumiLabel()
plot.SetNormFactor(True)
plot.Logy(True)
plot.SetTDRStyle()
plot.InitLegend()

plotr = root.HistogramDrawer()
plotr.SetLumi(lumi)
plotr.SetRatioStyle()
plotr.AddCMSLabel()
plotr.AddLumiLabel()
plotr.SetGrid()
#plotr.InitLegend(.15,.6,.5,.8)

plot.AddHistogram(hfake,'LF MC',root.kGjets)
plot.AddHistogram(hdata,'Data-(HF MC)',root.kData)
plot.Draw(args.outdir,'csvweight_shapes')

plotr.AddHistogram(hratio,'Data/MC',root.kExtra1,root.kBlack)
plotr.AddAdditional(hratioerr,'e2')
plotr.Draw(args.outdir,'csvweight_ratio')


fout = root.TFile(args.outdir+'/csvweight_fake.root','RECREATE')
fout.WriteTObject(hdata)
fout.WriteTObject(hfake)
fout.WriteTObject(hratio)
fout.Close()
