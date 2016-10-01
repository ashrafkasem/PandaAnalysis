#!/bin/bash

WD=$SUBMIT_WORKDIR

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
  tar -chzf 8011.tgz src python biglib bin lib objs test external # h = --dereference symlinks
  mv 8011.tgz $WD
fi

wget -O $WD/local.cfg $PANDA_CFG

cd ${CMSSW_BASE}/src/PandaAnalysis/SubMIT/inputs/
sed "s?XXXX?${SUBMIT_OUTDIR}?g" skim_tmpl.py > skim.py
cp skim.py $WD
chmod 775 ${WD}/skim.py

voms-proxy-init -voms cms
cp /tmp/x509up_u$UID $WD/x509up

cp ${CMSSW_BASE}/src/PandaAnalysis/SubMIT/inputs/exec.sh $WD

# input files for submission: 8011.tgz, skim.py, x509up, local.cfg. exec.sh is the executable
