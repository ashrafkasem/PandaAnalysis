from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
}

metFilter='metFilter==1'


presel = 'nJet>=1 && jet1Pt>100 && jet1IsTight && fabs(jet1Eta)<2.4'
cuts = {
    # analysis regions
    'signal'            : tAND(metFilter,tAND(presel,'pfmet>200 && dphipfmet>0.5 && (nLooseMuon+nLooseElectron+nLoosePhoton+nTau)==0 && fabs(calomet-pfmet)/pfmet<0.5')), 
    #'singlemuon'        : tAND(metFilter,tAND(presel,'pfUWmag>200 && dphipfUW>0.5 && (nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==1 && looseLep1IsTight==1 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160')),
    'singlemuon'        : tAND(metFilter,tAND(presel,'pfUWmag>200 && dphipfUW>0.5 && (nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==1 && looseLep1IsTight==1 && fabs(calomet-pfmet)/pfUWmag<0.5')),
    #'singleelectron'    : tAND(metFilter,tAND(presel,'pfUWmag>200 && dphipfUW>0.5 && (nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==1 && looseLep1IsTight==1 && looseLep1IsHLTSafe==1 && pfmet>50 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160')),
    #'singleelectron'    : tAND(metFilter,tAND(presel,'pfUWmag>200 && dphipfUW>0.5 && nLoosePhoton==0 && nTau==0 && nLooseLep==1 && looseLep1IsTight==1 && looseLep1IsHLTSafe==1 && abs(looseLep1PdgId)==11 && pfmet>50 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160')),
    'singleelectron'    : tAND(metFilter,tAND(presel,'pfUWmag>200 && dphipfUW>0.5 && nLoosePhoton==0 && nTau==0 && nLooseLep==1 && nTightLep>0 && abs(looseLep1PdgId)==11 && pfmet>50 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160')),
    'dimuon'            : tAND(metFilter,tAND(presel,'pfUZmag>200 && dphipfUZ>0.5 && (nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==2 && looseLep1IsTight==1 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5')),
    'dielectron'        : tAND(metFilter,tAND(presel,'pfUZmag>200 && dphipfUZ>0.5 && (nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==2 && looseLep1IsTight==1 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5')),
    'photon'            : tAND(metFilter,tAND(presel,'pfUAmag>200 && dphipfUA>0.5 && (nLooseMuon+nLooseElectron+nTau)==0 && nLoosePhoton==1 && loosePho1IsTight==1 && fabs(loosePho1Eta)<1.4442 && fabs(calomet-pfmet)/pfUAmag<0.5')),
}

weights = {
    # analysis weights
  'signal'         : '%f*normalizedWeight*sf_lep*sf_lepReco*sf_ewkV*sf_qcdV*sf_metTrig',
  'top'            : '%f*normalizedWeight*sf_lep*sf_lepReco*sf_ewkV*sf_qcdV',
  'w'              : '%f*normalizedWeight*sf_lep*sf_lepReco*sf_ewkV*sf_qcdV',
  'notag'          : '%f*normalizedWeight*sf_lep*sf_lepReco*sf_ewkV*sf_qcdV',
  'singleelectron' : '%f*normalizedWeight*sf_lep*sf_lepReco*sf_ewkV*sf_qcdV*sf_eleTrig',
  'singlemuon'     : '%f*normalizedWeight*sf_lep*sf_lepReco*sf_ewkV*sf_qcdV*sf_metTrig',
  'dielectron'     : '%f*normalizedWeight*sf_lep*sf_lepReco*sf_ewkV*sf_qcdV*sf_eleTrig',
  'dimuon'         : '%f*normalizedWeight*sf_lep*sf_lepReco*sf_ewkV*sf_qcdV*sf_metTrig',
  'photon'         : '%f*normalizedWeight*sf_ewkV*sf_qcdV*sf_pho*sf_phoTrig',
}

for x in ['singlemuontop','singleelectrontop']:
  weights[x] = weights['top']
for x in ['singlemuonw','singleelectronw']:
	weights[x] = weights['w']
