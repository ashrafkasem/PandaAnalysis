#!/bin/bash

export PANDA="${CMSSW_BASE}/src/PandaAnalysis"
export PANDA_PROD="${EOS2}/pandaprod/v5/" # can have multiple paths, separated by : 
export PANDA_CFG="http://snarayan.web.cern.ch/snarayan/eoscatalog/20161005.cfg"
export PANDA_FLATDIR="${HOME}/home000/panda/v4/"

#export SUBMIT_CFG="test"
export SUBMIT_CFG="prod"
export SUBMIT_LOGDIR="/work/sidn/panda/logs/"
export SUBMIT_WORKDIR="/work/sidn/panda/submit/"
export SUBMIT_OUTDIR="/mnt/hadoop/cms/store/user/snarayan/panda/v2/batch/"
export SUBMIT_TMPL="skim_tmpl.py"
