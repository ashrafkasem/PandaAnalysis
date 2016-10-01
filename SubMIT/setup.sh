#!/bin/bash

export PANDA_PROD="${EOS2}/pandaprod/v3/" # can have multiple paths, separated by : 
export PANDA_CFG="http://snarayan.web.cern.ch/snarayan/eoscatalog/20161001.cfg"

#export SUBMIT_CFG="test"
export SUBMIT_CFG="prod"
export SUBMIT_LOGDIR="/work/sidn/panda/logs/"
export SUBMIT_WORKDIR="/work/sidn/panda/submit/"
export SUBMIT_OUTDIR="/mnt/hadoop/cms/store/user/snarayan/panda/v0/batch/"
