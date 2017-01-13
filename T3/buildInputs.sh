#!/bin/bash

WD=$SUBMIT_WORKDIR
rm -rf $WD/*

doTar=0
while getopts ":t" opt; do
  case $opt in
    t)
      doTar=1
      ;;
  esac 
done

cd $CMSSW_BASE/
if [[ $doTar == 1 ]]; then
  echo "Tarring up CMSSW..."
  tar -chzf cmssw.tgz src python biglib bin lib objs test external # h = --dereference symlinks
  mv cmssw.tgz $WD
fi

wget -O $WD/local.cfg $PANDA_CFG
cp $WD/local.cfg $WD/local_all.cfg # useful for cataloging

cd ${CMSSW_BASE}/src/PandaAnalysis/T3/inputs/
sed "s?XXXX?${SUBMIT_OUTDIR}?g" ${SUBMIT_TMPL} > skim.py
cp skim.py $WD
chmod 775 ${WD}/skim.py

voms-proxy-init -voms cms
cp /tmp/x509up_u$UID $WD/x509up

cp ${CMSSW_BASE}/src/PandaAnalysis/T3/inputs/exec.sh $WD

cp -r $WD $SUBMIT_OUTDIR/workdir/

# input files for submission: cmssw.tgz, skim.py, x509up, local.cfg. exec.sh is the executable
