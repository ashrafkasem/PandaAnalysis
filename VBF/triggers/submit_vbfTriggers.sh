#!/bin/bash
function splitPath {
  local IFS=/
  local pieces
  # Disable glob expansion
  set -f
  pieces=( $@ ) 
  set +f
  #printf '%d\n' "${#pieces[@]}"
  #printf '%s\n' "${pieces[@]}"
  echo ${pieces[${#pieces[@]}-1]}
}

# name of the job
jobName=$1
# text file with list of file locations
listOfDataFiles=$2
# (electrons | muons)
flavor=$3

startDir=`pwd`
receptacle=/data/t3home000/${USER}/vbfTriggerPlots

mkdir -p ${receptacle}/${jobName}
echo "Making tarball of CMSSW: ${receptacle}/CMSSW_${jobName}.tgz"
cd $CMSSW_BASE
rm ${receptacle}/*.tgz
tar  --exclude-vcs -chzf "${receptacle}/CMSSW_${jobName}.tgz" src/PandaAnalysis src/PandaCore src/PandaTree  python biglib bin lib objs test external
to_run="/home/dhsu/bin/condor-run ./batch_vbfTriggers.sh"
libs="${receptacle}/CMSSW_${jobName}.tgz $CMSSW_BASE/src/PandaAnalysis/VBF/triggers/batch_vbfTriggers.sh"
jobNumber=0
cd $CMSSW_BASE/src
while read line
do
    cd ${receptacle}/${jobName}
    options=( $line )
    batchFile=${options[0]}
    batchFileBase=`splitPath ${batchFile}`
    outputFile="vbf_batchTree_triggerEff_${jobNumber}.root"
    echo "Setting up job with input files: ${libs}; output files ${outputFile}"
    args="--auxiliary-input ${libs} --output ${outputFile}"
    args="${args} --task-name ${jobName}"
    #args="${args} --job-args \"${flavor} ${batchFileBase} ${jobNumber} ${jobName} ${CMSSW_VERSION}\" "
    args="${args} --job-args \"${flavor} ${batchFile} ${jobNumber} ${jobName} ${CMSSW_VERSION}\" "
    #echo "$to_run $args"
    eval "$to_run $args"
    ((jobNumber++))
done < "${listOfDataFiles}"
cd ${startDir}
