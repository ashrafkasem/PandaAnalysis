#!/bin/bash

export PANDA="${CMSSW_BASE}/src/PandaAnalysis"
#export PANDA_CFG="http://t3serv001.mit.edu/~snarayan/histcatalog/20170205.cfg"
export PANDA_CFG="http://t3serv001.mit.edu/~snarayan/histcatalog/20170302_zptt.cfg"
export PANDA_FLATDIR="${HOME}/home000/store/panda/v_8026_0_0_scimitar/"
mkdir -p $PANDA_FLATDIR

#export SUBMIT_TMPL="skim_csv_tmpl.py"
export SUBMIT_TMPL="skim_scimitar_tmpl.py"
export SUBMIT_NAME="v_8026_0_0_scimitar"
export SUBMIT_WORKDIR="/data/t3serv014/snarayan/condor/"${SUBMIT_NAME}"/work/"
export SUBMIT_LOGDIR="/data/t3serv014/snarayan/condor/"${SUBMIT_NAME}"/logs/"
export SUBMIT_OUTDIR="/data/t3serv014/snarayan/panda/"${SUBMIT_NAME}"/batch/"
#export SUBMIT_OUTDIR="/mnt/hadoop/scratch/snarayan/panda/"${SUBMIT_NAME}"/batch/"
mkdir -p $SUBMIT_WORKDIR $SUBMIT_OUTDIR/locks/ $SUBMIT_LOGDIR

#private production 
export PRIVATE_LOGDIR="${HOME}/cms/logs/monotop_private_panda/"
export PRIVATE_PRODDIR="${HOME}/cms/hist/monotop_private_panda/"
export PRIVATE_CFGDIR="${HOME}/cms/condor/monotop_private_panda/"

# fitting
export PANDA_FIT=/data/t3serv014/snarayan/CMSSW_7_4_7/
export PANDA_FITSCAN=${PANDA_FLATDIR}/scans/
export PANDA_LIMITS=${PANDA_FLATDIR}/limits/
export PANDA_XSECS=/home/snarayan/cms/cmssw/analysis/MonoTop_Xsec/
