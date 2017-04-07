#!/usr/bin/env python

from sys import argv 
from os import getenv
import argparse
from collections import namedtuple
import numpy as np

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--sample',metavar='sample',type=str)
args = parser.parse_args()
argv = []
basedir = getenv('PANDA_FLATDIR')

import root_numpy as rnp
import ROOT as root
from PandaCore.Tools.Load import Load

Load('CanvasDrawer')
plot = root.CanvasDrawer()
plot.SetTDRStyle()
c = root.TCanvas()


if args.sample=='QCD':
	weightname = 'ptweight_qcd'
	fin = root.TFile.Open(basedir+'/QCD_evt8.root')
	cut = None
else:
	weightname = 'ptweight_binned'
	fin = root.TFile.Open(basedir+'/ZpTT.root')
	cut = 'fj1IsMatched==1 && fj1GenSize<1.44'
#fin = root.TFile.Open('/home/snarayan/home000/store/panda/v_8024_4_3//Diboson.root')
#weightname = 'normalizedWeight'
t = fin.Get('events')

NBINSX=10
NBINSY=40
DrawCfg = namedtuple('DrawCfg',['branch','label','lo','hi'])
def draw(var1,var2,histname):
    h2 = root.TH2F(histname,'',NBINSX,var1.lo,var1.hi,NBINSY,var2.lo,var2.hi)
    h2.GetXaxis().SetTitle(var1.label)
    h2.GetYaxis().SetTitle(var2.label)
    h2.GetYaxis().SetTitleOffset(1.5)
    xarr = rnp.tree2array(tree=t, branches=[var1.branch,var2.branch,weightname], selection=cut)
    arr = np.array([xarr[var1.branch],xarr[var2.branch]])
    arr = arr.transpose()
    rnp.fill_hist(hist=h2,array=arr,weights=xarr[weightname])
    for iX in xrange(1,NBINSX+1):
        h2.SetBinContent(iX,1,h2.GetBinContent(iX,0)
                              +h2.GetBinContent(iX,1))
        h2.SetBinContent(iX,NBINSY,h2.GetBinContent(iX,NBINSY+1)
                                   +h2.GetBinContent(iX,NBINSY))
    h2.SetLineColor(15)
    h2.SetFillColor(15)
    h2.SetMarkerStyle(20)
    return h2

pt = DrawCfg('fj1Pt','fatjet p_{T} [GeV]',250,1000)
mSD = DrawCfg('fj1MSD','fatjet m_{SD} [GeV]',0,500)
nSDConst = DrawCfg('fj1NSDConst','N_{SD constituents}',0,2000)
sdEFrac = DrawCfg('fj1SDEFrac100','Energy frac of 100 hardest SD particles',0.97,1.001)
nConst = DrawCfg('fj1NConst','N_{particles}',0,2000)
EFrac = DrawCfg('fj1EFrac100','Energy frac of 100 hardest particles',0.97,1.001)

for cfg1,cfg2 in [(pt,sdEFrac),(pt,nSDConst),(pt,EFrac),(pt,nConst),(pt,mSD)]:
#for cfg1,cfg2 in [(pt,sdEFrac)]:
    histname = '%s_%s'%(cfg1.branch,cfg2.branch)
    hist = draw(cfg1,cfg2,histname)
    c.Clear()
    #hist.Draw('COLZ TEXT')
    hist.Draw('VIOLIN')
    plot.AddCMSLabel(.16,.94)
    plot.SetLumi(36.6); plot.AddLumiLabel(True)
    plot.SetCanvas(c)
    plot.Draw(args.outdir+'/',args.sample+'_'+histname)
