from PandaCore.Tools.Misc import *

fixedmjj = "mjj"

cuts = {}
weights = {}
triggers = {}

eventsel = 'metfilter==1 && filterbadChCandidate==1 && filterbadPFMuon==1 && fabs(minJetMetDPhi_withendcap) > 0.5 && (fabs(caloMet-trueMet)/met) < 0.5 && n_tau==0 && n_bjetsMedium==0 && met>200'
#noid = tAND(eventsel,'jot1Eta*jot2Eta < 0 && jot1Pt>80. && jot2Pt>40. && fabs(jot1Eta)<4.7 && fabs(jot2Eta)<4.7') 
noid = tAND(eventsel,'jot1Eta*jot2Eta < 0 && jot1Pt>30. && jot2Pt>30. && fabs(jot1Eta)<4.7 && fabs(jot2Eta)<4.7 && (fabs(jot1Eta)<3||fabs(jot1Eta)>3.2)') 
baseline = noid

#vbf cuts
cuts['baseline'] = baseline
cuts['noid' ] = noid
cuts['dEtaCut'] = tAND(baseline,'fabs(jjDEta)>3.')
cuts['mjjCut'] = tAND(baseline, '%s>500'%fixedmjj)
cuts['dEtaAndMjjCut'] =tAND(baseline, 'fabs(jjDEta)>3 && %s>500'%fixedmjj)

#regions
cuts['signal'] = tAND(baseline,'n_loosepho==0 && n_looselep==0')
cuts['wmn'] = tAND(baseline,'n_loosepho==0 && n_looselep==1 && abs(lep1PdgId)==13 && n_tightlep>0 && mt<160')
cuts['wen'] = tAND(baseline,'n_loosepho==0 && n_looselep==1 && abs(lep1PdgId)==11 && n_tightlep>0 && mt<160 && trueMet>50')
cuts['zmm'] = tAND(baseline,'n_loosepho==0 && n_looselep == 2 && lep2PdgId*lep1PdgId==-169 && n_tightlep>0 && fabs(dilep_m-91)<30')
cuts['zee'] = tAND(baseline,'n_loosepho==0 && n_looselep == 2 && lep2PdgId*lep1PdgId==-121 && n_tightlep>0 && fabs(dilep_m-91)<30')
cuts['pho'] = tAND(baseline,'photonPt>175 && fabs(photonEta)<1.4442 && n_mediumpho==1 && n_loosepho==1 && n_looselep == 0')
cuts['zll'] = tOR(cuts['zmm'],cuts['zee'])
cuts['wlv'] = tOR(cuts['wmn'],cuts['wen'])

weights['signal'] = 'normalizedWeight*lepton_SF1*lepton_SF2*METTrigger*puWeight*topPtReweighting'
weights['zmm'] = tTIMES(weights['signal'],'tracking_SF1*tracking_SF2')
weights['wmn'] = tTIMES(weights['signal'],'tracking_SF1')
weights['zee'] = tTIMES(weights['signal'],'gsfTracking_SF1*gsfTracking_SF2')
weights['wen'] = tTIMES(weights['signal'],'gsfTracking_SF1')
weights['pho'] = 'normalizedWeight*PhoTrigger*topPtReweighting*photon_SF*puWeight'

triggers['met'] = "(triggerFired[10]==1 || triggerFired[11] == 1 || triggerFired[12] || triggerFired[13] == 1)" 
triggers['ele'] = "(triggerFired[0] || triggerFired[1] || triggerFired[2] || triggerFired[3] || triggerFired[4] || triggerFired[5] || triggerFired[26])"
triggers['pho'] = "(triggerFired[18] || triggerFired[19] || triggerFired[17] || triggerFired[5] || triggerFired[15] || triggerFired[16])"
