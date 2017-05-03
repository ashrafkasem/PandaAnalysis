#!/usr/bin/env python

from os import system,environ,readlink,path
from sys import exit,stdout,argv
from glob import glob
from PandaCore.Tools.condor import classad,htcondor

user = environ['USER']
fittingdir = environ['PANDA_FITTING']
scansdir = fittingdir+'/scans2/'
scramdir = environ['PANDA_FIT']
flatdir  = environ['PANDA_FLATDIR']

argv = []
import ROOT as root

sigfiles = []
nper = 5
force = False
#nper = 1

# sigfiles += glob(flatdir+'/ST*root')
# coupling_scan = False

#sigfiles += glob(flatdir+'/Scalar*root')
#coupling_scan = False

sigfiles += glob(flatdir+'/Vector*root')
coupling_scan = True

if coupling_scan:
    fweights = open('../signal_weights.dat')
    weights = [x.strip().replace('_nlo','').replace('rw_','') for x in fweights]
else:
    weights = ['nominal']


cfgheader = '%s %s'%(fittingdir,scansdir)
cfglines = []
mask = {}

for ff in sigfiles:
    f = ff.split('/')[-1].replace('.root','')
    model = ''
    mParams = ''
    mV = ''
    mChi = ''
    couplings = ['']
    if 'Vector' in f:
        replacements = {
            'Vector_MonoTop_NLO_Mphi-':'',
            '_gSM-0p25_gDM-1p0_13TeV-madgraph':'',
            '_Mchi-':'_',
            }
        for k,v in replacements.iteritems():
            f = f.replace(k,v)
        mV,mChi = map(int,f.split('_'))
        model = 'vector'
        for w in weights:
            cfg = '%s %i_%i %s'%(model,mV,mChi,w) 
            if force:
                mask[cfg] = True
            else:
                if len(glob(scansdir+'/'+model+'/'+w+'/higgsCombine*%i_%i*root'%(mV,mChi)))==0:
                    mask[cfg] = True
                else:
                    checkfilepath = glob(scansdir+'/'+model+'/'+w+'/higgsCombine*%i_%i*root'%(mV,mChi))[0]
                    f = root.TFile(checkfilepath)
                    t = f.Get('limit')
                    if not(t) or t.GetEntriesFast()<6:
                        print 'corrupt:',checkfilepath
                        mask[cfg] = True
                    else:
                        mask[cfg] = False
            cfglines.append(cfg)

exit(1)

to_submit = []
for i in xrange(len(cfglines)):
    cfg = cfglines[i]
    if mask[cfg]:
        to_submit.append(i)


fout = open('submit.cfg','w')
fout.write(cfgheader+'\n')
for l in cfglines:
    fout.write(l+'\n')
fout.close()

args_template = '%s --template correlated_tmpl.txt --cfg %s '%(scramdir,path.abspath('./submit.cfg'))

coll = htcondor.Collector()
schedd = htcondor.Schedd(coll.locate(htcondor.DaemonTypes.Schedd, 't3home000.mit.edu'))
base_job_properties = {
    "Cmd" : "runLimit2.sh",
    "WhenToTransferOutput" : "ON_EXIT",
    "ShouldTransferFiles" : "YES",
    "Requirements" : classad.ExprTree('UidDomain == "mit.edu" && Arch == "X86_64" && OpSysAndVer == "SL6"'),
    "AcctGroup" : "group_t3mit.urgent",
    "AccountingGroup" : "group_t3mit.urgent.snarayan",
    "OnExitHold" : classad.ExprTree("( ExitBySignal == true ) || ( ExitCode != 0 )"),
}
cluster_ad = classad.ClassAd()
for k,v in base_job_properties.iteritems():
    cluster_ad[k] = v
procs = []
for counter in xrange(len(to_submit)/nper+1):
    proc_ad = classad.ClassAd()
    proc_ad['UserLog'] = fittingdir + '/logs/%i.log'%counter
    proc_ad['Out'] = fittingdir + '/logs/%i.out'%counter
    proc_ad['Err'] = fittingdir + '/logs/%i.err'%counter
    args = args_template
    for i in to_submit[nper*counter:min(nper*(counter+1),len(to_submit))]:
        args += ' %i'%(i)
    proc_ad['Arguments'] = args
    procs.append((proc_ad,1))
print 'Submitting %i jobs'%(len(procs))
schedd.submitMany(cluster_ad,procs)

