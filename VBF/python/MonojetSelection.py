from PandaCore.Tools.Misc import *

fixedmjj = "mjj"

cuts = {}
weights = {}
triggers = {}

eventsel = 'metfilter==1 && filterbadChCandidate==1 && filterbadPFMuon==1 && fabs(minJetMetDPhi_withendcap)>0.5 && (fabs(caloMet-trueMet)/met)<0.5 && n_tau==0 && met>200'
noid = tAND(eventsel,'jet1Pt>100. && leadingJet_outaccp==0 && fabs(jet1Eta)<2.4') 
baseline = tAND(noid,'jet1isMonoJetIdNew==1')

#regions
cuts['signal'] = tAND(baseline,'n_loosepho==0 && n_looselep==0')
cuts['wmn'] = tAND(baseline,'n_loosepho==0 && n_looselep==1 && abs(lep1PdgId)==13 && n_tightlep>0 && mt<160')
#cuts['wmn'] = 'jet1isMonoJetIdNew==1&&leadingJet_outaccp==0&&1==1&&abs(jet1Eta)<2.5&&n_tau==0&&n_bjetsMedium==0&&n_looselep == 1 && abs(lep1PdgId)==13 && mt<160&&n_loosepho==0&&metfilter==1 && filterbadChCandidate==1 && filterbadPFMuon==1&&jet1Pt>100. &&(abs(caloMet-trueMet)/met) < 0.5&&n_tightlep > 0&&abs(minJetMetDPhi_withendcap) > 0.5&&met>200.0'
cuts['wen'] = tAND(baseline,'n_loosepho==0 && n_looselep==1 && abs(lep1PdgId)==11 && n_tightlep>0 && trueMet>50 && mt<160 && lep1IsHLTSafe==1')
cuts['zmm'] = tAND(baseline,'n_loosepho==0 && n_looselep == 2 && lep2PdgId*lep1PdgId==-169 && n_tightlep>0 && fabs(dilep_m-91)<30')
cuts['zee'] = tAND(baseline,'n_loosepho==0 && n_looselep == 2 && lep2PdgId*lep1PdgId==-121 && n_tightlep==1 && fabs(dilep_m-91)<30')
cuts['pho'] = tAND(baseline,'photonPt>175 && fabs(photonEta)<1.4442 && n_mediumpho==1 && n_loosepho==1 && n_looselep == 0')
cuts['zll'] = tOR(cuts['zmm'],cuts['zee'])
cuts['wlv'] = tOR(cuts['wmn'],cuts['wen'])

weights['signal'] = 'normalizedWeight*topPtReweighting*METTrigger'
weights['zmm'] = 'normalizedWeight*topPtReweighting*lepton_SF1*lepton_SF2*tracking_SF1*tracking_SF2*METTrigger'
weights['zee'] = 'normalizedWeight*topPtReweighting*lepton_SF1*lepton_SF2*gsfTracking_SF1*gsfTracking_SF2'
weights['wmn'] = 'normalizedWeight*topPtReweighting*lepton_SF1*tracking_SF1*METTrigger'
weights['wen'] = 'normalizedWeight*topPtReweighting*lepton_SF1*gsfTracking_SF1*EleTrigger'
weights['pho'] = 'normalizedWeight*PhoTrigger*topPtReweighting*photon_SF*PhoTrigger'

triggers['met'] = ' || '.join(['triggerFired[%i]'%x for x in [54,55,56,57,58,59,60,61,62,63,64,65]])
triggers['ele'] = ' || '.join(['triggerFired[%i]'%x for x in [7,8,11,12,13,42,43,44,45]])
triggers['pho'] = ' || '.join(['triggerFired[%i]'%x for x in [73,74,75,76,77]])
