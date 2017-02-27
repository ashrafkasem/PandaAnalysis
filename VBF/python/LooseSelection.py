from PandaCore.Tools.Misc import *

fixedmjj = "mjj"

cuts = {}
weights = {}
triggers = {}

eventsel = 'metfilter==1 && filterbadChCandidate==1 && filterbadPFMuon==1 && fabs(minJetMetDPhi_withendcap) > 0.5 && (fabs(caloMet-trueMet)/met) < 0.5 && n_tau==0 && n_bjetsMedium==0 && met>200'
#eventsel = 'metfilter==1 && filterbadChCandidate==1 && filterbadPFMuon==1 && fabs(minJetMetDPhi_withendcap) > 0.5 && n_tau==0 && n_bjetsMedium==0 && met>200'
#noid = tAND(eventsel,'jot1Eta*jot2Eta < 0 && jot1Pt>80. && jot2Pt>40. && fabs(jot1Eta)<4.7 && fabs(jot2Eta)<4.7') 
noid = tAND(eventsel,'jot1Eta*jot2Eta < 0 && jot1Pt>80. && jot2Pt>40. && fabs(jot1Eta)<4.7 && fabs(jot2Eta)<4.7 && (fabs(jot1Eta)<3||fabs(jot1Eta)>3.2)') 
baseline = noid

#regions
cuts['signal'] = tAND(baseline,'n_loosepho==0 && n_looselep==0')
cuts['wmn'] = tAND(baseline,'n_loosepho==0 && n_looselep==1 && abs(lep1PdgId)==13 && n_tightlep>0 && mt<160')
cuts['wen'] = tAND(baseline,'n_loosepho==0 && n_looselep==1 && abs(lep1PdgId)==11 && n_tightlep>0 && mt<160 && trueMet>50')
cuts['zmm'] = tAND(baseline,'n_loosepho==0 && n_looselep == 2 && lep2PdgId*lep1PdgId==-169 && n_tightlep>0 && fabs(dilep_m-91)<30')
cuts['zee'] = tAND(baseline,'n_loosepho==0 && n_looselep == 2 && lep2PdgId*lep1PdgId==-121 && n_tightlep>0 && fabs(dilep_m-91)<30')
cuts['pho'] = tAND(baseline,'photonPt>175 && fabs(photonEta)<1.4442 && n_mediumpho==1 && n_loosepho==1 && n_looselep == 0')
cuts['zll'] = tOR(cuts['zmm'],cuts['zee'])
cuts['wlv'] = tOR(cuts['wmn'],cuts['wen'])

weights['signal'] = '%f*normalizedWeight*METTrigger*topPtReweighting'
weights['zmm'] = '%f*normalizedWeight*lepton_SF1_id*lepton_SF1_iso*lepton_SF2_id*lepton_SF2_iso*METTrigger*topPtReweighting' 
weights['wmn'] = '%f*normalizedWeight*lepton_SF1_id*lepton_SF1_iso*METTrigger*topPtReweighting' 
weights['zee'] = '%f*normalizedWeight*lepton_SF1_v2*gsfTracking_SF1_v2*lepton_SF2_v2*gsfTracking_SF2_v2*topPtReweighting' 
weights['wen'] = '%f*normalizedWeight*lepton_SF1_v2*gsfTracking_SF1_v2*EleTrigger*topPtReweighting' 

triggers['met'] = ' || '.join(['triggerFired[%i]'%x for x in range(56,68)])
triggers['ele'] = ' || '.join(['triggerFired[%i]'%x for x in [7,8,11,12,13,43,44,45,46]])
triggers['pho'] = ' || '.join(['triggerFired[%i]'%x for x in range(75,80)]) 
