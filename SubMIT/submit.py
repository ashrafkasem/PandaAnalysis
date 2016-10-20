#!/usr/bin/env python

from sys import argv
from os import system,getenv,getuid

cmssw_base=getenv('CMSSW_BASE')
logpath=getenv('SUBMIT_LOGDIR')
workpath=getenv('SUBMIT_WORKDIR')
uid=getuid()

cfgpath = workpath+'/local.cfg'
cfgfile = open(cfgpath)
njobs = len(list(cfgfile))
nper = 5
njobs = njobs/nper+1

classad='''
universe = vanilla
executable = {0}/exec.sh
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = {0}/8011.tgz,{0}/local.cfg,{0}/skim.py,{0}/x509up
transfer_output_files = ""
input = /dev/null
output = {1}/$(Cluster)_$(Process).out
error = {1}/$(Cluster)_$(Process).err
log = {1}/$(Cluster)_$(Process).log
requirements = ( \ 
                  (OSGVO_OS_STRING == "RHEL 6" && HAS_CVMFS_cms_cern_ch) || \
                  GLIDEIN_REQUIRED_OS == "rhel6" || \
                  (GLIDEIN_Site == "MIT_CampusFactory" && (BOSCOGroup == "paus" || BOSCOGroup == "bosco_cms") && HAS_CVMFS_cms_cern_ch) \
                ) && \
                ( isUndefined(GLIDEIN_Entry_Name) || \
                  !stringListMember(GLIDEIN_Entry_Name,"CMS_T2_US_Nebraska_Red_op,CMS_T2_US_Nebraska_Red_gw1_op,CMS_T2_US_Nebraska_Red_gw2_op,CMS_T3_MX_Cinvestav_proton_work,CMS_T3_US_Omaha_tusker,CMSHTPC_T3_US_Omaha_tusker,Glow_US_Syracuse_condor,Glow_US_Syracuse_condor-ce01,Gluex_US_NUMEP_grid1,HCC_US_BNL_gk01,HCC_US_BNL_gk02,HCC_US_BU_atlas-net2,OSG_US_FIU_HPCOSGCE,OSG_US_Hyak_osg,OSG_US_UConn_gluskap,OSG_US_SMU_mfosgce",",") \ 
                ) && \ 
                ( isUndefined(GLIDEIN_Site) || !stringListMember(GLIDEIN_Site,"SU-OG,HOSTED_BOSCO_CE",",") )
rank = Mips
arguments = $(Process) {4}
use_x509userproxy = True
x509userproxy = /tmp/x509up_u{2}
#on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
+AccountingGroup = "analysis.sidn"
+AcctGroup = "analysis"
+ProjectName = "CpDarkMatterSimulation"
+DESIRED_Sites = "T3_US_SDSC,T2_FR_GRIF_LLR,T3_US_Omaha,T2_CH_CERN_AI,T2_IT_Bari,T2_CH_CERN,T2_CH_CSCS,T2_UA_KIPT,T2_IN_TIFR,T2_FR_IPHC,T2_IT_Rome,T2_UK_London_Brunel,T2_EE_Estonia,T2_US_Florida,T2_US_Wisconsin,T2_HU_Budapest,T2_DE_RWTH,T2_BR_UERJ,T2_ES_IFCA,T2_DE_DESY,T2_US_Caltech,T2_TW_Taiwan,T0_CH_CERN,T1_RU_JINR_Disk,T2_UK_London_IC,T2_US_Nebraska,T2_ES_CIEMAT,T3_US_Princeton,T2_PK_NCP,T2_CH_CERN_T0,T3_US_FSU,T3_KR_UOS,T3_IT_Perugia,T3_US_Minnesota,T2_TR_METU,T2_AT_Vienna,T2_US_Purdue,T3_US_Rice,T3_HR_IRB,T2_BE_UCL,T3_US_FIT,T2_UK_SGrid_Bristol,T2_PT_NCG_Lisbon,T1_ES_PIC,T3_US_JHU,T2_IT_Legnaro,T2_RU_INR,T3_US_FIU,T3_EU_Parrot,T2_RU_JINR,T2_IT_Pisa,T3_UK_ScotGrid_GLA,T3_US_MIT,T2_CH_CERN_HLT,T2_MY_UPM_BIRUNI,T1_FR_CCIN2P3,T2_FR_GRIF_IRFU,T2_FR_CCIN2P3,T2_PL_Warsaw,T3_AS_Parrot,T2_US_MIT,T2_BE_IIHE,T2_RU_ITEP,T1_CH_CERN,T3_CH_PSI,T3_IT_Bologna"
queue {3}
'''.format(workpath,logpath,uid,njobs,nper)


with open(logpath+'/condor.jdl','w') as jdlfile:
  jdlfile.write(classad)

system('condor_submit %s/condor.jdl'%logpath)
