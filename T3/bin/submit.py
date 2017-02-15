#!/usr/bin/env python

from sys import argv
from os import system,getenv,getuid
from time import time

cmssw_base=getenv('CMSSW_BASE')
logpath=getenv('SUBMIT_LOGDIR')
workpath=getenv('SUBMIT_WORKDIR')
uid=getuid()

cfgpath = workpath+'/local.cfg'
cfgfile = open(cfgpath)
njobs = len(list(cfgfile))
nper = 10
#nper = 50
njobs = njobs/nper + 1
#njobs = 1

now = int(time())
frozen_cfgpath = cfgpath.replace('local','local_%i'%now)
system('cp %s %s'%(cfgpath,frozen_cfgpath)) # freeze the config file

classad='''
universe = vanilla
executable = {0}/exec.sh
should_transfer_files = YES
#when_to_transfer_output = ON_EXIT
transfer_input_files = {0}/cmssw.tgz,{5},{0}/skim.py,{0}/x509up
#transfer_output_files = dummy.txt
input = /dev/null
output = {1}/$(Cluster)_$(Process).out
error = {1}/$(Cluster)_$(Process).err
log = {1}/$(Cluster)_$(Process).log
requirements = UidDomain == "mit.edu" && Arch == "X86_64" && OpSysAndVer == "SL6"
arguments = "$(Process) {4}"
use_x509userproxy = True
x509userproxy = /tmp/x509up_u{2}
accounting_group = group_cmsuser.snarayan
queue {3}
'''.format(workpath,logpath,uid,njobs,nper,frozen_cfgpath)

print classad

with open(logpath+'/condor.jdl','w') as jdlfile:
  jdlfile.write(classad)

system('condor_submit %s/condor.jdl'%logpath)
