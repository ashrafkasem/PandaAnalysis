#!/usr/bin/env python 

import argparse
from sys import argv,exit
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--process',metavar='process',type=str)
parser.add_argument('--selection',metavar='selection',type=str,default='monojet')
basedir = '/mnt/hadoop/scratch/snarayan/kfactors/'
args = parser.parse_args()

sname = argv[0]
argv = []


from os import getenv
import ROOT as root
from math import sqrt
from array import array
from PandaCore.Tools.Load import Load
from PandaCore.Tools.root_interface import draw_hist, read_files
from PandaCore.Tools.Misc import *
import PandaAnalysis.VBF.PandaSelection as sel

Load('PandaCoreTools')
Load('PandaCoreDrawers')

# create some global variables

plotr = root.HistogramDrawer()
plotr.SetRatioStyle()
plotr.SetGrid()
plotr.InitLegend(.7,.65,.88,.9)
plotr.SetAutoRange(False)

f_lo = basedir+args.process+'_lo.root'
f_nlo = basedir+args.process+'_nlo.root'

labels = {
        'a' : '#gamma',
        'z' : 'Z',
        'w' : 'W',
        'monofatjet' : 'fatjet p_{T}>100 GeV',
        'inclusive' : 'Inclusive',
        'monojet' : 'jet p_{T}>100 GeV',
        'monojethigh' : 'jet p_{T}>200 GeV',
        'mass' : 'fat jet p_{T}>100 GeV, m_{SD}>50',
        'tag' : 'fat jet p_{T}>100 GeV, m_{SD}>110, #tau_{32}<0.7',
        'tagonly' : 'fat jet p_{T}>100 GeV, #tau_{32}<0.7',
        }

proc_cuts = {
        'a' : '1==1',
        'z' : 'v_m>60 && v_m<120',
        'w' : '1==1',
        }

nlo_corrections = {
        'a' : '1',
        'z' : '3 / 2', # 3=flavors, 2=BR(vv)/BR(ll)
        'w' : '3',
        }

selections = {
        'inclusive' : '1==1',
        'monojet' : 'jpt_1>100',
        'monojethigh' : 'jpt_1>200',
        'monofatjet' : 'fjpt>100',
        'mass' : 'fjmsd>50 && fjpt>100',
        'tag'  : 'fjmsd>110 && fjt3t2<0.7 && fjpt>100',
        'tagonly'  : 'fjt3t2<0.7 && fjpt>100',
        }

bins = array('f',[175,225,275,325,375,425,500,600,700,1000])
nbins = len(bins)-1
hbase = root.TH1D('dummy','',nbins,bins)
hbase.GetXaxis().SetTitle('%s p_{T} [GeV]'%(labels[args.process]))
hbase.GetYaxis().SetTitle('k-factor')

ptstr = 'med_pt' if args.process=='a' else 'v_pt'

def draw(order):
    if order=='nlo':
        f = f_nlo
        weightstr = tTIMES('normalizedWeight',nlo_corrections[args.process])
    else:
        f = f_lo
        weightstr = 'normalizedWeight'
    cutstr = tAND(selections[args.selection],proc_cuts[args.process])

    h = hbase.Clone('h_'+order)
    xarr = read_files([f],[ptstr,weightstr],cutstr,treename='Events')
    draw_hist(h,xarr,[ptstr],weightstr)
    print order,h.Integral(),cutstr

    return h

def build_ratio(hnum,hden):
    hratio = hnum.Clone()
    hratio.Divide(hden)
    return hratio


h_lo = draw('lo')
h_nlo = draw('nlo')

h_kfac = build_ratio(h_nlo,h_lo)
h_kfac.SetMaximum(2.5)
h_kfac.SetMinimum(0.8)
h_err = h_kfac.Clone()
h_kfac.SetLineColor(root.kRed)
h_kfac.SetLineWidth(3)
h_err.SetLineColor(root.kGray)
h_err.SetFillColor(root.kGray)

plotr.cd()

plotr.AddCMSLabel()
plotr.AddPlotLabel('#splitline{%s NLO/LO}{%s}'%(labels[args.process],labels[args.selection])
                   ,.18,.7,False,42,.06)
plotr.AddHistogram(h_err,'stat. unc.',root.kExtra1,root.kRed,'e2')
plotr.AddHistogram(h_kfac,'NLO/LO',root.kExtra2,root.kRed,'hist')

plotname = 'kfac_%s_%s'%(args.process,args.selection)
plotr.Draw(args.outdir+'/',plotname)
