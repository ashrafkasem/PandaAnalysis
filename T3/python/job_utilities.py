#!/usr/bin/env python

from re import sub
from sys import exit
from os import system,getenv,path
from time import clock,time
import json

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import *
import PandaCore.Tools.job_management as cb

sname = 'T3.job_utilities'
data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'

_stopwatch = time() 
def print_time(label):
    global _stopwatch
    now_ = time()
    PDebug(sname+'.print_time:'+str(time()),
           '%.3f s elapsed performing "%s"'%((now_-_stopwatch),label))
    _stopwatch = now_

def copy_local(long_name):
    replacements = {
                r'\${EOS}':'root://eoscms.cern.ch//store/user/snarayan',
                r'\${EOS2}':'root://eoscms.cern.ch//store/group/phys_exotica',
                r'\${CERNBOX}':'root://eosuser//eos/user/s/snarayan',
                r'\${CERNBOXB}':'root://eosuser//eos/user/b/bmaier',
            }
    full_path = long_name
    for k,v in replacements.iteritems():
        full_path = sub(k,v,full_path)
    PInfo(sname,full_path)

    panda_id = long_name.split('/')[-1].split('_')[-1].replace('.root','')
    input_name = 'input_%s.root'%panda_id

    # if the file is cached locally, why not use it?
    local_path = full_path.replace('root://xrootd.cmsaf.mit.edu/','/mnt/hadoop/cms')
    if path.isfile(local_path):
        # apparently SmartCached files can be corrupted...
        ftest = root.TFile(local_path)
        if ftest and not(ftest.IsZombie()):
            PInfo(sname+'.copy_local','Opting to read locally')
            return local_path

    # rely on pxrdcp for local and remote copies
    # default behavior: drop PF candidates
    cmd = "pxrdcp %s %s '!pfCandidates'"%(full_path,input_name)
    PInfo(sname+'.copy_local',cmd)

    system(cmd)
            
    if path.isfile(input_name):
        PInfo(sname+'.copy_local','Successfully copied to %s'%(input_name))
        return input_name
    else:
        PError(sname+'.copy_local','Failed to copy %s'%input_name)
        return None


def cleanup(fname):
    ret = system('rm -f %s'%(fname))
    if ret:
        PError(sname+'.cleanup','Removal of %s exited with code %i'%(fname,ret))
    else:
        PInfo(sname+'.cleanup','Removed '+fname)
    return ret


def hadd(good_inputs):
    good_outputs = ' '.join([x.replace('input','output') for x in good_inputs])
    cmd = 'hadd -f output.root ' + good_outputs
    ret = system(cmd)    
    if not ret:
        PInfo(sname+'.hadd','Merging exited with code %i'%ret)
    else:
        PError(sname+'.hadd','Merging exited with code %i'%ret)


def drop_branches(to_drop=None, to_keep=None):
    # remove any irrelevant branches from the final tree.
    # this MUST be the last step before stageout or you 
    # run the risk of breaking something
    
    if not to_drop and not to_keep:
        return 0

    if to_drop and to_keep:
        PError(sname+'.drop_branches','Can only provide to_drop OR to_keep')
        return 0

    f = root.TFile('output.root','UPDATE')
    t = f.FindObjectAny('events')
    n_entries = t.GetEntriesFast() # to check the file wasn't corrupted
    if to_drop:
        if type(to_drop)==str:
            t.SetBranchStatus(to_drop,False)
        else:
            for b in to_drop:
                t.SetBranchStatus(b,False)
    elif to_keep:
        t.SetBranchStatus('*',False)
        if type(to_keep)==str:
            t.SetBranchStatus(to_keep,True)
        else:
            for b in to_keep:
                t.SetBranchStatus(b,True)
    t_clone = t.CloneTree()
    f.WriteTObject(t_clone,'events','overwrite')
    f.Close()

    # check that the write went okay
    f = root.TFile('output.root')
    if f.IsZombie():
        PError(sname+'.drop_branches','Corrupted file trying to drop '+to_drop)
        return 1 
    t_clone = f.FindObjectAny('events')
    if (n_entries==t_clone.GetEntriesFast()):
        return 0
    else:
        PError(sname+'.drop_branches','Corrupted tree trying to drop '+to_drop)
        return 2


def stageout(outdir,outfilename):
    mvargs = 'mv $PWD/output.root %s/%s'%(outdir,outfilename)
    PInfo(sname,mvargs)
    ret = system(mvargs)
    system('rm *.root')
    if not ret:
        PInfo(sname+'.stageout','Move exited with code %i'%ret)
    else:
        PError(sname+'.stageout','Move exited with code %i'%ret)
        return ret
    if not path.isfile('%s/%s'%(outdir,outfilename)):
        PError(sname+'.stageout','Output file is missing!')
        ret = 1
    return ret


def write_lock(outdir,outfilename,processed):
    flock = open(outdir+'/locks/'+outfilename.replace('.root','.lock'),'w')
    for k,v in processed.iteritems():
        flock.write(v+'\n')
    flock.close()


def classify_sample(full_path, isData):
    if not isData:
        if any([x in full_path for x in ['Vector_','Scalar_']]):
            return root.PandaAnalyzer.kSignal
        elif any([x in full_path for x in ['ST_','ZprimeToTT']]):
            return root.PandaAnalyzer.kTop
        elif 'EWKZ2Jets' in full_path:
            return root.PandaAnalyzer.kZEWK
        elif 'EWKW' in full_path:
            return root.PandaAnalyzer.kWEWK
        elif 'ZJets' in full_path or 'DY' in full_path:
            return root.PandaAnalyzer.kZ
        elif 'WJets' in full_path:
            return root.PandaAnalyzer.kW
        elif 'GJets' in full_path:
            return root.PandaAnalyzer.kA
        elif 'TTJets' in full_path or 'TT_' in full_path:
            return root.PandaAnalyzer.kTT
    return root.PandaAnalyzer.kNone


def add_json(skimmer, json_path):
    with open(json_path) as jsonFile:
        payload = json.load(jsonFile)
        for run_str,lumis in payload.iteritems():
            run = int(run_str)
            for l in lumis:
                skimmer.AddGoodLumiRange(run,l[0],l[1])

def run_PandaAnalyzer(skimmer, isData, input_name):
    # some common stuff that doesn't need to be configured
    # read the inputs
    try:
        fin = root.TFile.Open(input_name)
        tree = fin.FindObjectAny("events")
        weight_table = fin.FindObjectAny('weights')
        hweights = fin.FindObjectAny("hSumW")
    except:
        PError(sname+'.run_PandaAnalyzer','Could not read %s'%input_name)
        return False # file open error => xrootd?
    if not tree:
        PError(sname+'.run_PandaAnalyzer','Could not recover tree in %s'%input_name)
        return False
    if not hweights:
        PError(sname+'.run_PandaAnalyzer','Could not recover hweights in %s'%input_name)
        return False
    if not weight_table:
        weight_table = None

    output_name = input_name.replace('input','output')
    skimmer.SetDataDir(data_dir)
    if isData:
        add_json(skimmer, data_dir+'/certs/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt')

    rinit = skimmer.Init(tree,hweights,weight_table)
    if rinit:
        PError(sname+'.run_PandaAnalyzer','Failed to initialize %s!'%(input_name))
        return False 
    skimmer.SetOutputFile(output_name)

    # run and save output
    skimmer.Run()
    skimmer.Terminate()

    ret = path.isfile(output_name)
    if ret:
        PInfo(sname+'.run_PandaAnalyzer','Successfully created %s'%(output_name))
        return True
    else:
        PError(sname+'.run_PandaAnalyzer','Failed in creating %s!'%(output_name))
        return False


def main(to_run, processed, fn):
    print_time('loading')
    for f in to_run.files:
        input_name = copy_local(f)
        print_time('copy %s'%input_name)
        if input_name:
            success = fn(input_name,(to_run.dtype!='MC'),f)
            print_time('analyze %s'%input_name)
            if success:
                processed[input_name] = f
            if input_name[:5] == 'input': # if this is a local copy
                cleanup(input_name)
            print_time('remove %s'%input_name)
    
    if len(processed)==0:
        PWarning(sname+'.main', 'No successful outputs!')
        exit(1)


