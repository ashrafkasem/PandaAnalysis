#!/usr/bin/env python 

import argparse
from sys import argv
from os import getenv
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--process',metavar='process',type=str,default=None)
parser.add_argument('--process1',metavar='process1',type=str,default=None)
parser.add_argument('--process2',metavar='process2',type=str,default=None)
parser.add_argument('--region1',metavar='region1',type=str)
parser.add_argument('--region2',metavar='region2',type=str)
args = parser.parse_args()
basedir = getenv('PANDA_FLATDIR')+'/'
if args.process:
    args.process1 = args.process    
    args.process2 = args.process    

sname = argv[0]
argv = []
import ROOT as root
from math import sqrt
from collections import namedtuple 
from array import array
from PandaCore.Tools.Load import Load
from PandaCore.Tools.root_interface import draw_hist, read_tree
from PandaCore.Tools.Misc import *
import PandaAnalysis.VBF.PandaSelection as sel

Load('PandaCoreTools')
Load('PandaCoreDrawers')

# create some global variables
plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.InitLegend()
plot.SetLumi(35.8)

plotr = root.HistogramDrawer()
plotr.SetRatioStyle()
plotr.InitLegend(.7,.65,.88,.9)
plotr.SetLumi(35.8)
plotr.SetAutoRange(False)

repl = {}

repl['up'] = {
    'pfmet' : 'pfmetUp',
    'jot1Pt' : 'jot1PtUp',
    'jot2Pt' : 'jot2PtUp',
    'jot1Eta' : 'jot1EtaUp',
    'jot2Eta' : 'jot2EtaUp',
    'jot12Mass' : 'jot12MassUp',
    'jot12DEta' : 'jot12DEtaUp',
    'jot12DPhi' : 'jot12DPhiUp',
}
repl['down'] = {
    'pfmet' : 'pfmetDown',
    'jot1Pt' : 'jot1PtDown',
    'jot2Pt' : 'jot2PtDown',
    'jot1Eta' : 'jot1EtaDown',
    'jot2Eta' : 'jot2EtaDown',
    'jot12Mass' : 'jot12MassDown',
    'jot12DEta' : 'jot12DEtaDown',
    'jot12DPhi' : 'jot12DPhiDown',
}

labels = {
        'WJets' : 'W',
        'ZJets' : 'Z',
        'ZtoNuNu' : 'Z',
        }

f_in1 = root.TFile.Open(basedir+args.process1+'.root')
t_in1 = f_in1.Get('events')
f_in2 = root.TFile.Open(basedir+args.process2+'.root')
t_in2 = f_in2.Get('events')

hbase = root.TH1D('dummy','',1,0,2)
hbase.GetXaxis().SetTitle('1')
hbase.GetYaxis().SetTitle('Events')

cut1 = {'nominal' : sel.cuts[args.region1]}
cut2 = {'nominal' : sel.cuts[args.region2]}
weight1 = {'nominal' : sel.weights[args.region1]%35800}
weight2 = {'nominal' : sel.weights[args.region2]%35800}

for s in ['up','down']:
    repl_ = repl[s]
    for d in [cut1,cut2,weight1,weight2]:
        d[s] = d['nominal']
        for k,v in repl_.iteritems():
            d[s] = d[s].replace(k,v)
        d[s] = d[s].replace('dphipfmetUp','dphipfmet')
        d[s] = d[s].replace('dphipfmetDown','dphipfmet')

def draw(s):
    h1 = hbase.Clone(s+'_1')
    h2 = hbase.Clone(s+'_2')
    xarr1 = read_tree(t_in1,['1',weight1[s]],cut1[s])
    xarr2 = read_tree(t_in2,['1',weight2[s]],cut2[s]) 
    draw_hist(h1,xarr1,['1'],weight1[s])
    draw_hist(h2,xarr2,['1'],weight2[s])
    h1.SetMaximum(h1.GetMaximum()*2)
    h2.SetMaximum(h2.GetMaximum()*2)
    return {1:h1,2:h2}

def build_ratio(hnum,hden):
    hratio = hnum.Clone()
    hratio.Divide(hden)
    return hratio

nominal = draw('nominal')
up = draw('up')
down = draw('down')

for h in [nominal,up,down]:
    for i in [1,2]:
        PInfo(sname,'%30s => %f'%(str(h[i]), h[i].Integral()))

label1 = ''
if args.process1 in labels:
    label1 += labels[args.process1]
else:
    label1 += args.process1
if args.region1 == 'signal':
    label1 += ' (SR)'
elif 'electron' in args.region1:
    label1 += ' (e'
    if 'di' in args.region1:
        label1 += 'e)'
    else:
        label1 += ')'
elif 'muon' in args.region1:
    label1 += ' (#mu'
    if 'di' in args.region1:
        label1 += '#mu)'
    else:
        label1 += ')'
label2 = ''
if args.process2 in labels:
    label2 += labels[args.process2]
else:
    label2 += args.process2
if args.region2 == 'signal':
    label2 += ' (SR)'
elif 'electron' in args.region2:
    label2 += ' (e'
    if 'di' in args.region2:
        label2 += 'e)'
    else:
        label2 += ')'
elif 'muon' in args.region2:
    label2 += ' (#mu'
    if 'di' in args.region2:
        label2 += '#mu)'
    else:
        label2 += ')'

# main plotting functions
def plot_ratio():
   
    h_nominal = build_ratio(nominal[1],nominal[2]) 
    h_up = build_ratio(up[1],up[2]) 
    h_up = build_ratio(h_up,h_nominal)
    h_down = build_ratio(down[1],down[2]) 
    h_down = build_ratio(h_down,h_nominal)
    h_nominal = build_ratio(h_nominal,h_nominal)

    h_up.Divide(h_nominal)
    h_down.Divide(h_nominal)
    for ib in xrange(1,h_nominal.GetNbinsX()+1):
        val = h_nominal.GetBinContent(ib)
        if val==0: continue
        h_nominal.SetBinContent(ib,1)
        h_nominal.SetBinError(ib,h_nominal.GetBinError(ib)/val)
    for h in [h_nominal,h_up,h_down]:
        h.SetMaximum(1.5)
        h.SetMinimum(0.8)

    h_nominal.GetYaxis().SetTitle('Uncertainty')

    plotr.cd()
    plotr.Reset()

    plotr.AddCMSLabel()
    plotr.AddLumiLabel(True)
    plotr.AddPlotLabel('%s / %s'%(label1,label2),.18,.77,False,42,.06)
    plotr.AddHistogram(h_nominal,'Nominal (stat.)',root.kData)
    plotr.AddHistogram(h_up,'JEC Up',root.kExtra2,root.kRed,'hist')
    plotr.AddHistogram(h_down,'JEC Down',root.kExtra3,root.kBlue,'hist')

    plotname = 'jec_ratio_%s_%s_%s'%(args.process1,args.region1,args.region2)
    plotr.Draw(args.outdir+'/',plotname)

def plot_yield(region):
    plot.cd()
    plot.Reset()
    plot.AddCMSLabel()
    plot.AddLumiLabel(True)
    plot.AddHistogram(nominal[region],'Nominal (stat.)',root.kData)
    plot.AddHistogram(up[region],'JEC Up',root.kExtra2,root.kRed,'hist')
    plot.AddHistogram(down[region],'JEC Down',root.kExtra3,root.kBlue,'hist')

    process = args.process1 if region==1 else args.process2
    label_ = label1 if region==1 else label2
    plot.AddPlotLabel(label_,.18,.77,False,42,.05)
    plotname = 'jec_yieldo_%s_%s'%(process,args.region1 if region==1 else args.region2)
    plot.Draw(args.outdir+'/',plotname)

plot_ratio()
plot_yield(1)
plot_yield(2)
