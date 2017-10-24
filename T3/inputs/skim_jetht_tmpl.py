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
import PandaAnalysis.T3.job_utilities as utils

Load('PandaAnalyzer')
data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'



def fn(input_name, isData, full_path):
    
    PInfo(sname+'.fn','Starting to process '+input_name)
    # now we instantiate and configure the analyzer
    skimmer = root.PandaAnalyzer()
    skimmer.isData=isData
    skimmer.SetFlag('fatjet',False)
    skimmer.SetFlag('genOnly',True)
    skimmer.SetPreselectionBit(root.PandaAnalyzer.kRecoil50)
    processType = utils.classify_sample(full_path, isData)
    skimmer.processType=processType 

    return utils.run_PandaAnalyzer(skimmer, isData, input_name)


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


if __name__ == "__main__":
    sample_list = cb.read_sample_config('local.cfg',as_dict=False)
    to_run = None #sample_list[which]
    for s in sample_list:
        if which==s.get_id():
            to_run = s
            break
    if not to_run:
        PError(sname,'Could not find a job for PROCID=%i'%(which))
        exit(3)

    outdir = 'XXXX' # will be replaced when building the job
    outfilename = to_run.name+'_%i.root'%(submit_id)
    processed = {}
    
    utils.main(to_run, processed, fn)

    utils.hadd(list(processed))
    utils.print_time('hadd')

    ret = utils.stageout(outdir,outfilename)
    utils.print_time('stageout')
    if not ret:
        utils.write_lock(outdir,outfilename,processed)
        utils.print_time('create lock')
    else:
        exit(-1*ret)

    exit(0)
