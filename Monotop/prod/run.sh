#!/bin/bash

label=$1
cfgName=$2
outpath=$3
scramdir=$4/src/

echo $label $cfgName $outpath

pwd
rm *root

#executable=prod.py
executable=mcCondor.py

cd ${scramdir}/
eval `scramv1 runtime -sh`
cd -

echo RUNNING ON $HOSTNAME
voms-proxy-init -voms cms
#cp ${scramdir}/MitPanda/Monotop/prod/x509up_u67051 .
#export X509_USER_PROXY=${PWD}/x509up_u67051
#echo "X509_USER_PROXY = $X509_USER_PROXY"

pandadir=${scramdir}/PandaProd/Producer/cfg/

cat $cfgName

for f in $(cat $cfgName); do
  echo "copying $f"
  xrdcp ${f} .
done 

echo "##############################################################"

touch local.cfg
for f in $(ls *.root); do
  echo "file:${PWD}/${f}" >> local.cfg
done

cat local.cfg

echo "##############################################################"

ls

echo "##############################################################"

cmsRun ${pandadir}/${executable} filelist=local.cfg outfile=${PWD}/pandatree_${label}.root

echo "##############################################################"

mv ${PWD}/panda_${label}.root ${outpath}

rm -r *root local.cfg
