#!/bin/bash

export PATH=${PATH}:${CMSSW_BASE}/src/PandaCore/bin/

export PANDA="${CMSSW_BASE}/src/PandaAnalysis"
export PANDA_CFG="http://t3serv001.mit.edu/~bmaier/stuff/all_008.cfg" 
export PANDA_FLATDIR="/mnt/hadoop/scratch/bmaier/panda/008_v1/flat/"
mkdir -p $PANDA_FLATDIR

#export SUBMIT_TMPL="skim_noegm_tmpl.py"
export SUBMIT_TMPL="skim_wlnhbb_tmpl.py"
#export SUBMIT_TMPL="skim_scimitar_tmpl.py"
export SUBMIT_NAME="v_008_v1"
export SUBMIT_WORKDIR="/data/t3home000/bmaier/panda/"${SUBMIT_NAME}"/work/"
export SUBMIT_LOGDIR="/data/t3home000/bmaier/panda/panda/"${SUBMIT_NAME}"/logs/"
export SUBMIT_OUTDIR="/mnt/hadoop/scratch/bmaier/panda/"${SUBMIT_NAME}"/batch/"
mkdir -p $SUBMIT_WORKDIR $SUBMIT_OUTDIR/locks/ $SUBMIT_LOGDIR
