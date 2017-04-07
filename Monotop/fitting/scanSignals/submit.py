#!/usr/bin/env python

from os import system,environ
from sys import exit,stdout
from glob import glob

user = environ['USER']
fittingdir = environ['PANDA_FITTING']
scansdir = fittingdir+'/scans/'
scramdir = environ['PANDA_FIT']
flatdir  = environ['PANDA_FLATDIR']

sigfiles = []
#sigfiles += glob(flatdir+'/Scalar*root')
sigfiles += glob(flatdir+'/Vector*root')

coupling_scan = True

to_submit = {}

for ff in sigfiles:
    f = ff.split('/')[-1].replace('.root','')
    model = ''
    mParams = ''
    mV = ''
    mChi = ''
    couplings = ['']
    couplings = []
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
        if coupling_scan:
            couplings += [x.split('/')[-1] for x in glob(scansdir+'g*')]
    else:
        replacements = {
               'Scalar_MonoTop_LO_Mphi-':'',
               '_13TeV-madgraph':'',
               '_Mchi-':'_',
            }
        for k,v in replacements.iteritems():
            f = f.replace(k,v)
        mV,mChi = map(int,f.split('_'))
        model = '--isRes'
    mParams = '%i_%i'%(mV,mChi)
    for coupling in couplings:
        if coupling in to_submit:
            to_submit[coupling].append(mParams)
        else:
            to_submit[coupling] = [mParams]

if coupling_scan:
    for k in to_submit:
        to_submit[k] = '+'.join(to_submit[k])


for coupling,mParams in to_submit.iteritems():
    if coupling_scan:
        if coupling!='':
            logpath = fittingdir+'/logs/%s.'%(coupling)
        else:
            logpath = fittingdir+'/logs/nominal.'
    else:
        logpath = fittingdir+'/logs/%s_%s.'%(mParams,coupling)
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
                         Queue 1'''.format(logpath,fittingdir,scramdir,model,mParams,coupling,user)
    system('echo "%s" | condor_submit'%(condorJDLString))
