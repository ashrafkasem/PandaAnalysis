#!/usr/bin/env python

from os import system,environ
from sys import exit,stdout,argv
from time import sleep


user = environ['USER']
cmssw = environ['CMSSW_BASE']
cfgName=argv[1]
nPerJob = 50
xrd='root://xrootd.cmsaf.mit.edu/'

histdir = '/home/%s/cms/hist/monotop_private_pandatree/%s/'%(user,cfgName)
cfgdir = '/home/%s/cms/condor/monotop_private_panda/%s/'%(user,cfgName)
logdir = '/home/%s/cms/logs/monotop_private_panda/%s/'%(user,cfgName)
for dire in [histdir,cfgdir,logdir]:
  system('rm -rf '+dire)
  system('mkdir -p '+dire)

def submit(l):
  print "submitting",l,cfgName
  condorJDLString = '''
Executable                 = run.sh
Universe                   = vanilla
requirements = UidDomain == \\\"mit.edu\\\" && Arch == \\\"X86_64\\\" && OpSysAndVer == \\\"SL6\\\" && Machine != \\\"t3home000.mit.edu\\\"
Error                      = %s/\$(Process).err
Output                     = %s/\$(Process).out
Log                        = %s/\$(Process).log
Arguments                  = \$(Process) %s/\$(Process).cfg %s %s
should_transfer_files      = YES
when_to_transfer_output    = ON_EXIT
GetEnv                     = True
accounting_group           = group_cmsuser.%s
Queue %i'''%(logdir,logdir,logdir,cfgdir,histdir,cmssw,user,l)
  system('echo "%s" | condor_submit'%condorJDLString)


with open("cfg/%s.txt"%(cfgName)) as cfgFile:
  filesToProcess = list(cfgFile)


nJobs = len(filesToProcess)/nPerJob+1
for iJ in xrange(nJobs):
  filesThisJob = filesToProcess[nPerJob*iJ:nPerJob*(iJ+1)]
  with open(cfgdir+'%i.cfg'%(iJ),'w') as thisCfg:
    for f in filesThisJob:
      thisCfg.write(xrd+f)


submit(nJobs)
#submit(10)
