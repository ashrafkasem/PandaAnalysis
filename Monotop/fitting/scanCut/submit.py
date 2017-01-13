#!/usr/bin/env python

from os import system,environ
from sys import exit,stdout

#cuts = ['top_ecf_bdt>%.2f'%(0.1+0.05*x) for x in range(18)]
cuts = ['%.2f'%(0.1+0.05*x) for x in range(18)]
user = environ['USER']
outdir = environ['PANDA_FITSCAN']
scramdir = environ['PANDA_FIT']

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
