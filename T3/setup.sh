#!/bin/bash

export PATH=${PATH}:${CMSSW_BASE}/src/PandaCore/bin/

export PANDA="${CMSSW_BASE}/src/PandaAnalysis"
export PANDA_CFG="http://t3serv001.mit.edu/~snarayan/histcatalog/test_t2.cfg" 
export PANDA_FLATDIR="/home/snarayan/home000/store/panda/v_008_v0/"
mkdir -p $PANDA_FLATDIR

export SUBMIT_TMPL="skim_wlnhbb_tmpl.py"
#export SUBMIT_TMPL="skim_gghbb_tmpl.py"
export SUBMIT_NAME="v_008_v0"
export SUBMIT_WORKDIR="/scratch/snarayan/jobs/"${SUBMIT_NAME}"/work/"
export SUBMIT_LOGDIR="/scratch/snarayan/jobs/"${SUBMIT_NAME}"/logs/"
export SUBMIT_LOCKDIR="/scratch/snarayan/jobs/"${SUBMIT_NAME}"/locks/"
export SUBMIT_OUTDIR="/mnt/hadoop/scratch/snarayan/panda/"${SUBMIT_NAME}"/batch/"
mkdir -p $SUBMIT_WORKDIR $SUBMIT_OUTDIR $SUBMIT_LOGDIR $SUBMIT_LOCKDIR


# testing
export SUBMIT_CONFIG=T2
