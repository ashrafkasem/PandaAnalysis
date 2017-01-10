from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
}

metFilter='metFilter==1'
topTagSF = '%f*(fj1IsMatched==1)+%f*(fj1IsMatched==0)'%(1.,1.)
#ak4bTagSF = 'sf_btag0*(isojetNBtags==0)+sf_btag1*(isojetNBtags==1)+1*(isojetNBtags>1)'
ak4bTagSF = '0.95*sf_csvWeightB*sf_csvWeightM*sf_sjcsvWeightB*sf_sjcsvWeightM'
#ak4bTagSF = '1'

presel = 'nFatjet==1 && fj1Pt>250 && TMath::Abs(fj1Eta)<2.4 && 50<fj1MSD'
cuts = {
    'signal_nobtag'     : tAND(metFilter,tAND(presel,'pfmet>250 && dphipfmet>1.1 && (nLooseMuon+nLooseElectron+nLoosePhoton+nTau)==0')), # turning off btags
    'signal'            : tAND(metFilter,tAND(presel,'pfmet>250 && dphipfmet>1.1 && (nLooseMuon+nLooseElectron+nLoosePhoton+nTau)==0 && fj1MaxCSV>0.46 && isojetNBtags==0')),
    'singlemuontop'     : tAND(metFilter,tAND(presel,'pfUWmag>250 && (nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==1 && looseLep1IsTight==1 && fj1MaxCSV>0.46 && isojetNBtags==1')),
    'singleelectrontop' : tAND(metFilter,tAND(presel,'pfUWmag>250 && (nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==1 && looseLep1IsTight==1 && fj1MaxCSV>0.46 && isojetNBtags==1 && pfmet>40')),
    'singlemuonw'       : tAND(metFilter,tAND(presel,'pfUWmag>250 && (nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==1 && looseLep1IsTight==1 && fj1MaxCSV<0.46 && isojetNBtags==0')),
    'singleelectronw'   : tAND(metFilter,tAND(presel,'pfUWmag>250 && (nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==1 && looseLep1IsTight==1 && fj1MaxCSV<0.46 && isojetNBtags==0 && pfmet>40')),
    'singlemuon'        : tAND(metFilter,tAND(presel,'pfUWmag>250 && (nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==1 && looseLep1IsTight==1')),
    'singleelectron'    : tAND(metFilter,tAND(presel,'pfUWmag>250 && (nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==1 && looseLep1IsTight==1 && pfmet>40')),
    'dimuon'            : tAND(metFilter,tAND(presel,'pfUZmag>250 && (nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==2 && looseLep1IsTight==1 && 60<diLepMass && diLepMass<120')),
    'dielectron'        : tAND(metFilter,tAND(presel,'pfUZmag>250 && (nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==2 && looseLep1IsTight==1 && 60<diLepMass && diLepMass<120')),
    'photon'            : tAND(metFilter,tAND(presel,'pfUAmag>250 && (nLooseMuon+nLooseElectron+nTau)==0 && nLoosePhoton==1 && loosePho1IsTight==1')),
}

weights = {
  'signal_nobtag'    : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_ewkV*sf_qcdV*sf_tt8TeV',topTagSF),'1'),
  'signal'    : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_ewkV*sf_qcdV*sf_tt8TeV',topTagSF),ak4bTagSF),
  'top'       : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_ewkV*sf_qcdV*sf_tt8TeV',topTagSF),ak4bTagSF),
  'w'         : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_ewkV*sf_qcdV*sf_tt8TeV',topTagSF),ak4bTagSF),
  'notag'     : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_ewkV*sf_qcdV*sf_tt8TeV',topTagSF),'1'),
  'singleelectron' : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_ewkV*sf_qcdV*sf_tt8TeV',topTagSF),'1'),
  'singlemuon'     : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_ewkV*sf_qcdV*sf_tt8TeV',topTagSF),'1'),
}

def addTrig(region,weight):
  if 'electron' in region:
    return tTIMES('sf_eleTrig',weight)
  elif 'photon' in region:
    return tTIMES('sf_phoTrig',weight)
  else:
    return tTIMES('sf_metTrig',weight)

for x in ['singlemuontop','singleelectrontop']:
  weights[x] = addTrig(x,weights['top'])
for x in ['singlemuonw','singleelectronw']:
  weights[x] = addTrig(x,weights['w'])
for x in ['dimuon','dielectron']:
  #if x=='dielectron':
    weights[x] = addTrig(x,tTIMES('sf_lep',weights['notag']))
  #else:
  #  weights[x] = weights['notag']
for x in ['photon']:
  #weights[x] = addTrig(x,weights['notag'])
  weights[x] = tTIMES('sf_pho',weights['notag'])

for r in ['signal','top','w','singlemuontop','singleelectrontop','singlemuonw','singleelectronw']:
  for shift in ['BUp','BDown','MUp','MDown']:
    for cent in ['sf_btag','sf_sjbtag']:
      weights[r+'_'+cent+shift] = sub(cent+'0',cent+'0'+shift,sub(cent+'1',cent+'1'+shift,weights[r]))
