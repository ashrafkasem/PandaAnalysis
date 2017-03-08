#!/usr/bin/env python

from os import system,environ
from sys import exit,stdout
from glob import glob

user = environ['USER']
limitdir = environ['PANDA_LIMITS']
scramdir = environ['PANDA_FIT']
flatdir  = environ['PANDA_FLATDIR']

sigfiles = glob(flatdir+'/Vector*root')

iC=0
for ff in sigfiles:
  f = ff.split('/')[-1].replace('.root','')
  model = ''
  mV = ''
  mChi = ''
  if 'Vector' in f:
    replacements = {
      'Vector_MonoTop_NLO_Mphi-':'',
      '_gSM-0p25_gDM-1p0_13TeV-madgraph':'',
      '_Mchi-':'_',
      }
    for k,v in replacements.iteritems():
      f = f.replace(k,v)
    mV,mChi = map(int,f.split('_'))
    model = '--isFCNC'
  else:
    mV = f.split('_')[1].split('-')[1]
    mChi = '100'
    model = '--isRes'
  logpath = limitdir+'/logs/%i.'%(iC)
  condorJDLString = '''Executable = runLimit.sh
Universe  = vanilla
requirements = UidDomain == \\\"mit.edu\\\" && Arch == \\\"X86_64\\\" && OpSysAndVer == \\\"SL6\\\"
Error = {0}err
Output  = {0}out
Log = {0}log
Arguments = {1} {2} {3} {4} {5}
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
GetEnv = True
accounting_group = group_cmsuser.{6}
Queue 1'''.format(logpath,limitdir,scramdir,model,mV,mChi,user)
  system('echo "%s" | condor_submit'%(condorJDLString))
  iC+=1
