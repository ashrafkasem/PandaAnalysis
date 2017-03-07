#!/usr/bin/env python

from re import sub
from sys import argv,exit
from os import system,getenv,path
from time import clock,time
import json

which = int(argv[1])
submit_id = int(argv[2])
sname = argv[0]
argv=[]

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import *
import PandaCore.Tools.job_management as cb
import PandaAnalysis.Tagging.cfg_v8 as tagcfg

now = int(time())
Load('PandaAnalyzer')
data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'

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

    # xrdcp if remote, copy if local 
    # xrdcp seems faster than pxrdcp, even if pf candidates are dropped
    if 'root://' in full_path:
        system('xrdcp %s %s'%(full_path,input_name))
    else:
        system('cp %s %s'%(full_path,input_name))

    # rely on pxrdcp for local and remote copies
    # default behavior: drop PF candidates
    # system("pxrdcp %s %s '!pfCandidates'"%(full_path,input_name))

            
    if path.isfile(input_name):
        PInfo(sname+'.copy_local','Successfully copied to %s'%(input_name))
        return input_name
    else:
        PError(sname+'.copy_local','Failed to copy %s'%input_name)
        return None


def fn(input_name,isData,full_path):
    start=clock()
    
    PInfo(sname+'.fn','Starting to process '+input_name)
    # now we instantiate and configure the analyzer
    skimmer = root.PandaAnalyzer()
    skimmer.isData=isData
    skimmer.SetFlag('firstGen',True)
    skimmer.SetFlag('pfCands',True)
    #skimmer.SetPreselectionBit(root.PandaAnalyzer.kRecoil)
    skimmer.SetPreselectionBit(root.PandaAnalyzer.kFatjet)
    processType=root.PandaAnalyzer.kNone
    if not isData:
        if any([x in full_path for x in ['ST_','Vector_','ZprimeToTT']]):
            processType=root.PandaAnalyzer.kTop
        elif 'ZJets' in full_path or 'DY' in full_path:
            processType=root.PandaAnalyzer.kZ
        elif 'WJets' in full_path:
            processType=root.PandaAnalyzer.kW
        elif 'GJets' in full_path:
            processType=root.PandaAnalyzer.kA
        elif 'TTJets' in full_path or 'TT_' in full_path:
            processType=root.PandaAnalyzer.kTT
    skimmer.processType=processType 

    # read the inputs
    try:
        fin = root.TFile.Open(input_name)
        tree = fin.FindObjectAny("events")
        hweights = fin.FindObjectAny("hSumW")
    except:
        PError(sname+'.fn','Could not read %s'%input_name)
        return False # file open error => xrootd?
    if not tree:
        PError(sname+'.fn','Could not recover tree in %s'%input_name)
        return False
    if not hweights:
        PError(sname+'.fn','Could not recover hweights in %s'%input_name)
        return False

    output_name = input_name.replace('input','output')
    skimmer.SetDataDir(data_dir)
    if isData:
        with open(data_dir+'/certs/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt') as jsonFile:
            payload = json.load(jsonFile)
            for run_str,lumis in payload.iteritems():
                run = int(run_str)
                for l in lumis:
                    skimmer.AddGoodLumiRange(run,l[0],l[1])
    skimmer.SetOutputFile(output_name)
    skimmer.Init(tree,hweights)

    # run and save output
    skimmer.Run()
    skimmer.Terminate()

    ret = path.isfile(output_name)
    if ret:
        PInfo(sname+'.fn','Successfully created %s in %.2f sec'%(output_name,(clock()-start)/1000.))
        return True
    else:
        PError(sname+'.fn','Failed in creating %s!'%(output_name))
        return False


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


def add_bdt():
    # now run the BDT
    Load('TMVABranchAdder')
    ba = root.TMVABranchAdder()
    ba.treename = 'events'
    ba.defaultValue = -1.2
    ba.presel = 'fj1ECFN_2_4_20>0'
    for v in tagcfg.variables:
        ba.AddVariable(v[0],v[2])
    for v in tagcfg.formulae:
        ba.AddFormula(v[0],v[2])
    for s in tagcfg.spectators:
        ba.AddSpectator(s[0])
    ba.BookMVA('top_ecf_bdt',data_dir+'/trainings/top_ecfbdt_v8_BDT.weights.xml')
    ba.RunFile('output.root')


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


if __name__ == "__main__":
    sample_list = cb.read_sample_config('local.cfg',as_dict=False)
    to_run = sample_list[which]
    outdir = 'XXXX' # will be replaced when building the job
    outfilename = to_run.name+'_%i.root'%(submit_id)
    processed = {}

    for f in to_run.files:
        input_name = copy_local(f)
        if input_name:
            success = fn(input_name,(to_run.dtype!='MC'),f)
            if success:
                processed[input_name] = f
            cleanup(input_name)
    
    if len(processed)==0:
        exit(1)

    hadd(list(processed))
    add_bdt()
    if drop_branches(to_keep=['fj1*','top_ecf_bdt','*Number']):
        exit(2)
    ret = stageout(outdir,outfilename)
    if not ret:
        write_lock(outdir,outfilename,processed)
    else:
        exit(-1*ret)

    exit(0)
