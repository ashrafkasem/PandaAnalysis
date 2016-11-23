#!/usr/bin/env python

import argparse
from sys import argv

VERBOSE=False
parser = argparse.ArgumentParser(description='fit stuff')
parser.add_argument('--disc',metavar='disc',type=str,default='top_ecfv8_bdt')
parser.add_argument('--cut',metavar='cut',type=float,default=0.6)
parser.add_argument('--useMass',metavar='useMass',type=str,default='FALSE')
args = parser.parse_args()
argv=[]


from math import sqrt
import ROOT as root
from PandaCore.Tools.Misc import *

basedir = '/home/snarayan/home000/store/panda/v8/sf/'
nbins = 1

fdata = root.TFile(basedir+'photon__data.root'); tdata = fdata.Get('data')
fgjets = root.TFile(basedir+'photon__gjets.root'); tgjets = fgjets.Get('gjets')

c = root.TCanvas()

def getIntAndErr(t,cut,weight):
  hone = root.TH1F('h1','h1',1,-2,2)
  hone.Sumw2()
  t.Draw('1>>h1',tTIMES(cut,weight))
  val = hone.GetBinContent(1)
  err = hone.GetBinError(1)
  return val,err

def calcEffAndErr(p,perr,f,ferr):
  eff_ = p/(p+f)
  err_ = pow( perr * f / pow(p+f,2) , 2 )
  err_ += pow( ferr * p / pow(p+f,2) , 2 )
  err_ = sqrt(err_)
  print p,perr,f,ferr,eff_,err_
  return eff_,err_

def calcSFAndErr(d,derr,m,merr):
  sf_ = d/m
  err_ = pow(derr/m,2)
  err_ += pow(d*merr/pow(m,2),2)
  err_ = sqrt(err_)
  return sf_,err_

if args.useMass.upper()=='FALSE':
  if 'tau' in args.disc:
    passcut='%s<%f && mSD>50'
    failcut='%s>%f && mSD>50'
  else:
    passcut='%s>%f && mSD>50'
    failcut='%s<%f && mSD>50'
else:
  if 'tau' in args.disc:
    passcut='%s<%f && mSD>110 && mSD<210'
    failcut='!(%s<%f && mSD>110 && mSD<210) && mSD>50'
  else:
    passcut='%s>%f && mSD>110 && mSD<210'
    failcut='!(%s>%f && mSD>110 && mSD<210) && mSD>50'

# first get the yields for the tagging cut
datapass,datapasserr = getIntAndErr(tdata,passcut%(args.disc,args.cut),'1')
datafail,datafailerr = getIntAndErr(tdata,failcut%(args.disc,args.cut),'1')
gjetspass,gjetspasserr = getIntAndErr(tgjets,passcut%(args.disc,args.cut),'weight')
gjetsfail,gjetsfailerr = getIntAndErr(tgjets,failcut%(args.disc,args.cut),'weight')

if VERBOSE:
  print datapass,datafail,datapasserr,datafailerr
  print gjetspass,gjetsfail,gjetspasserr,gjetsfailerr

mceff,mcerr = calcEffAndErr(gjetspass,gjetspasserr,gjetsfail,gjetsfailerr)
dataeff,dataerr = calcEffAndErr(datapass,datapasserr,datafail,datafailerr)
tagsf, tagsferr = calcSFAndErr(dataeff,dataerr,mceff,mcerr)

print '%.3g +/- %.2g'%(tagsf,tagsferr)
