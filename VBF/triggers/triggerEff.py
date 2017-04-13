#!/usr/bin/env python

from sys import argv,exit
from os import getenv
from array import array
from math import sqrt
import argparse
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--cut',metavar='cut',type=str,default=None)
parser.add_argument('--label',metavar='label',type=str)
args = parser.parse_args()
basedir = getenv('PANDA_FLATDIR')+'/'

sname = argv[0]
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.root_interface import draw_hist, read_tree
from PandaCore.Tools.Misc import *
import PandaAnalysis.VBF.PandaSelection as sel

Load('PandaCoreTools')
Load('PandaCoreDrawers')

plot = root.HistogramDrawer()
plot.SetLumi(35.8)
plot.Ratio(True)
plot.FixRatio(0.1)
plot.SetTDRStyle()
plot.InitLegend(0.5,0.8,0.88,0.9,2)
plot.SetRatioLabel('W/Z')
plot.SetAutoRange(False)

w_base_cut = 'metFilter==1 && egmFilter==1 && jot1Eta*jot2Eta<0 && jot1Pt>80 && jot2Pt>40 && fabs(jot1Eta)<4.7 && fabs(jot2Eta)<4.7 && (fabs(jot1Eta)<3||fabs(jot1Eta)>3.2) && nTau==0 && nLoosePhoton==0 && nLooseLep==1 && looseLep1IsTight==1 && abs(looseLep1PdgId)==13 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160 && dphipfUW>0.5'
z_base_cut = 'metFilter==1 && egmFilter==1 && jot1Eta*jot2Eta<0 && jot1Pt>80 && jot2Pt>40 && fabs(jot1Eta)<4.7 && fabs(jot2Eta)<4.7 && (fabs(jot1Eta)<3||fabs(jot1Eta)>3.2) && nTau==0 && nLooseElectron==0 && nLoosePhoton==0 && nLooseMuon==2 && nTightLep>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5 && mT<160 && dphipfUZ>0.5'
trigger = '(trigger&1)!=0'

recoil_bins = array('f',[60,120,180,240,320,500,1000])
#recoil_bins = array('f',[60,90,120,150,180,210,240,280,320,400,500,750,1000])
hbase = root.TH1D('dummy','',len(recoil_bins)-1,recoil_bins)
hbase.GetXaxis().SetTitle('MET no #mu [GeV]')
hbase.GetYaxis().SetTitle('Efficiency')
hbase.SetMaximum(1.2)
hbase.SetMinimum(0.2)

counter=0

fIn = root.TFile(basedir+'/SingleMuon.root')
events = fIn.Get('events')

def get_hist(cut):
    h = hbase.Clone()
    xarr = read_tree(events,['pfmetnomu'],cut)
    draw_hist(h,xarr,['pfmetnomu'],None)
    return h

def run_eff(base_cut):
    cut = tAND(base_cut,args.cut)
    hden = get_hist(cut)
    hnum = get_hist(tAND(cut,trigger))
    hratio = hnum.Clone()
    for ib in xrange(1,len(recoil_bins)):
        vnum = hnum.GetBinContent(ib)
        enum = hnum.GetBinError(ib)
        vden = hden.GetBinContent(ib)
        if vden==0:
            vden = 1
        hratio.SetBinContent(ib,vnum/vden)
        hratio.SetBinError(ib,enum/vden)
    
    return hratio
    

h_w_eff = run_eff(w_base_cut)
h_z_eff = run_eff(z_base_cut)
h_one = hbase.Clone()
for ib in xrange(1,len(recoil_bins)): h_one.SetBinContent(ib,1)
h_one.SetLineColor(root.kGray)
h_one.SetLineStyle(2)

plot.AddHistogram(h_w_eff,'Single-#mu',root.kData,root.kBlue,'el')
plot.AddHistogram(h_z_eff,'Di-#mu',root.kExtra1,root.kRed,'el')
plot.AddAdditional(h_one,'hist')

plot.Draw(args.outdir,args.label)
