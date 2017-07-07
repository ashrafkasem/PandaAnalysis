pwd
mkdir ${5}
mv CMSSW_${4}.tgz ${5}
cd ${5}
tar -zxf "CMSSW_${4}.tgz"
cd src
pwd
ls  PandaAnalysis/VBF/triggers
cmsenv
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:.:PandaAnalysis/VBF/triggers
root -b -l -q PandaAnalysis/VBF/triggers/triggerEff.C+\(\"$1\",\"$2\",\"$3\"\)
mv vbf_batchTree_triggerEff_${3}.root ../..
