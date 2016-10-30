#!/usr/bin/env python

from os import system,environ
from sys import exit,stdout

deta = [1.,1.5,2.,2.5,3,3.25,3.5,3.75,4,4.25,4.5,4.75]
dphi = [3.14*x/10 for x in range(1,11)] 
mjj = range(0,1200,150)
pt1 = [80+20*x for x in range(5)]
pt2 = [40+10*x for x in range(9)]

cuts = []
#cuts += ['jjDEta>%.2f&&mjj>%i'%(e,m) for e in deta for m in mjj]
cuts += ['jot1Pt>%.2f&&jot2Pt>%.2f'%(x,y) for x in pt1 for y in pt2]
#cuts += ['jjDEta>%.2f&&jot1Pt>%.2f'%(x,y) for x in deta for y in pt1]
#cuts += ['jjDEta>%.2f&&jot2Pt>%.2f'%(x,y) for x in deta for y in pt2]
#cuts += ['mjj>%i&&jot1Pt>%.2f'%(x,y) for x in mjj for y in pt1]
#cuts += ['mjj>%i&&jot2Pt>%.2f'%(x,y) for x in mjj for y in pt2]
#cuts += ['jjDEta>%.2f&&jjDPhi<%.2f'%(e,p) for e in deta for p in dphi]
#cuts += ['jjDEta>%.2f&&fabs(minJetMetDPhi_withendcap)>%.2f'%(e,p) for e in deta for p in dphi]
user = environ['USER']
outdir = environ['PANDA_VBFSCAN']
scramdir = environ['PANDA_VBFFIT']

iC=0
for cut in cuts:
  print cut
  logpath = outdir+'/logs/%i.'%(iC)
  condorJDLString = '''Executable = scanCut.sh
Universe  = vanilla
requirements = UidDomain == \\\"mit.edu\\\" && Arch == \\\"X86_64\\\" && OpSysAndVer == \\\"SL6\\\"
Error = {0}err
Output  = {0}out
Log = {0}log
Arguments = '{1}' {2} {3}
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
GetEnv = True
accounting_group = group_cmsuser.{4}
Queue 1'''.format(logpath,cut,scramdir,outdir,user)
  with open('condor.jdl','w') as jdlFile:
    jdlFile.write(condorJDLString)
  system('echo "%s" | condor_submit'%(condorJDLString))
  iC+=1
