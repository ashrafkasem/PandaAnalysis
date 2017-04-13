from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
}

metFilter='metFilter==1 && egmFilter==1'


#presel = 'jot1Eta*jot2Eta<0 && jot1Pt>80 && jot2Pt>40 && fabs(jot1Eta)<4.7 && fabs(jot2Eta)<4.7 && (fabs(jot1Eta)<3||fabs(jot1Eta)>3.2) && nTau==0  && jot12Mass>1500 && fabs(jot12DEta)>4.2 && fabs(jot12DPhi)<1.3'
presel = 'jot1Eta*jot2Eta<0 && jot1Pt>80 && jot2Pt>40 && fabs(jot1Eta)<4.7 && fabs(jot2Eta)<4.7 && (fabs(jot1Eta)<3||fabs(jot1Eta)>3.2) && nTau==0'
cuts = {
    'signal'             : tAND(metFilter,tAND(presel,'pfmet>200   && nLooseLep==0 && nLoosePhoton==0 ')), 
    'singlemuon'         : tAND(metFilter,tAND(presel,'pfUWmag>200 && nLoosePhoton==0 && nLooseLep==1 && looseLep1IsTight==1 && abs(looseLep1PdgId)==13  && mT<160')),
    'singleelectron'     : tAND(metFilter,tAND(presel,'pfUWmag>200 && nLoosePhoton==0 && nLooseLep==1 && looseLep1IsTight==1 && abs(looseLep1PdgId)==11  && mT<160 && pfmet>50')),
    'dimuon'             : tAND(metFilter,tAND(presel,'pfUZmag>200 && nLooseElectron==0 && nLoosePhoton==0 && nLooseMuon==2 && nTightLep>0 && 60<diLepMass && diLepMass<120 ')),
    'dielectron'         : tAND(metFilter,tAND(presel,'pfUZmag>200 && nLooseMuon==0 && nLoosePhoton==0 && nLooseElectron==2 && nTightLep>0 && 60<diLepMass && diLepMass<120 ')),
}
'''
cuts = {
    'signal'             : tAND(metFilter,tAND(presel,'pfmet>200   && dphipfmet>0.5 && nLooseLep==0 && nLoosePhoton==0 && fabs(calomet-pfmet)/pfmet<0.5')), 
    'singlemuon'         : tAND(metFilter,tAND(presel,'pfUWmag>200 && dphipfUW>0.5 && nLoosePhoton==0 && nLooseLep==1 && looseLep1IsTight==1 && abs(looseLep1PdgId)==13 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160')),
    'singleelectron'     : tAND(metFilter,tAND(presel,'pfUWmag>200 && dphipfUW>0.5 && nLoosePhoton==0 && nLooseLep==1 && looseLep1IsTight==1 && abs(looseLep1PdgId)==11 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160 && pfmet>50')),
    'dimuon'             : tAND(metFilter,tAND(presel,'pfUZmag>200 && dphipfUZ>0.5 && nLooseElectron==0 && nLoosePhoton==0 && nLooseMuon==2 && nTightLep>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5')),
    'dielectron'         : tAND(metFilter,tAND(presel,'pfUZmag>200 && dphipfUZ>0.5 && nLooseMuon==0 && nLoosePhoton==0 && nLooseElectron==2 && nTightLep>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5')),
}
'''


weights = {
  'signal'         : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_metTrig',
  'w'              : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV',
  'z'              : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV',
  'photon'         : '%f*sf_pu*normalizedWeight*sf_ewkV*sf_qcdV*sf_pho*sf_phoTrig *sf_qcdV2j', # add the additional 2-jet kfactor
}

for x in ['dimuon','dielectron','singlemuon','singleelectron']:
	if 'electron' in x:
	  weights[x] = tTIMES(weights['z'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['z'],'sf_metTrig')
