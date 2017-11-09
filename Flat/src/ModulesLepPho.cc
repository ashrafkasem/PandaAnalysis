#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define EGMSCALE 1

using namespace panda;
using namespace std;

void PandaAnalyzer::SimpleLeptons() {
  looseLep1PdgId=-1, looseLep2PdgId=-1;
  //electrons
  for (auto& ele : event.electrons) {
    float pt = ele.pt(); float eta = ele.eta(); float aeta = fabs(eta);
    if (pt<10 || aeta>2.5) continue;
    if (!ele.veto) continue;
    if (!ElectronIP(ele.eta(),ele.dxy,ele.dz)) continue;
    ele.setPtEtaPhiM(pt,eta,ele.phi(),511e-6);
    unsigned iL=gt->nLooseElectron;
    bool isFake   = ele.hltsafe;
    bool isMedium = ele.medium;
    bool isTight  = ele.tight;
    if (isTight) gt->nTightElectron++;
    int eleSelBit            = kLoose;
    if (isFake  ) eleSelBit |= kFake;
    if (isMedium) eleSelBit |= kMedium;
    if (isTight ) eleSelBit |= kTight;
    gt->electronPt[iL]           = pt;
    gt->electronEta[iL]          = eta;
    gt->electronPhi[iL]          = ele.phi();
    gt->electronSfLoose[iL]      = GetCorr(cEleLoose, eta, pt);
    gt->electronSfTight[iL]      = GetCorr(cEleTight, eta, pt);
    gt->electronSfUnc[iL]        = GetError(cEleMedium, eta, pt);
    gt->electronSfReco[iL]       = GetCorr(cEleReco, eta, pt);
    gt->electronSelBit[iL]       = eleSelBit;
    gt->electronPdgId[iL]        = ele.charge*-11;
    gt->electronTripleCharge[iL] = ele.tripleCharge;
    looseLeps.push_back(&ele);
    matchLeps.push_back(&ele);
    matchEles.push_back(&ele);
    gt->nLooseElectron++;
    if (gt->nLooseElectron>=2) break;
  }
  // muons
  for (auto& mu : event.muons) {
    float pt = mu.pt(); float eta = mu.eta(); float aeta = fabs(eta);
    if (pt<10 || aeta>2.4) continue;
    if (!mu.loose) continue;
    mu.setPtEtaPhiM(pt,eta,mu.phi(),0.106);
    bool isFake   = mu.tight  && mu.combIso()/mu.pt() < 0.4 && mu.chIso/mu.pt() < 0.4;
    bool isMedium = mu.medium && mu.combIso()/mu.pt() < 0.15;
    bool isTight  = mu.tight  && mu.combIso()/mu.pt() < 0.15;
    if (isTight) gt->nTightMuon++;
    int muSelBit            = kLoose;
    if (isFake)   muSelBit |= kFake;
    if (isMedium) muSelBit |= kMedium;
    if (isTight)  muSelBit |= kTight;
    unsigned iL=gt->nLooseMuon;
    gt->muonPt[iL]                   = pt;
    gt->muonEta[iL]                  = eta;
    gt->muonPhi[iL]                  = mu.phi();
    gt->muonSfLoose[iL]              = GetCorr(cMuLooseID, TMath::Abs(mu.eta()), mu.pt()) * GetCorr(cMuLooseIso, TMath::Abs(mu.eta()), mu.pt());
    gt->muonSfTight[iL]              = GetCorr(cMuTightID, TMath::Abs(mu.eta()), mu.pt())  * GetCorr(cMuTightIso, TMath::Abs(mu.eta()), mu.pt());
    gt->muonSfUnc[iL]                = GetError(cMuMediumID , TMath::Abs(mu.eta()), mu.pt());
    gt->muonSfReco[iL]               = GetCorr(cMuReco, mu.eta());
    gt->muonSelBit[iL]               = muSelBit;
    gt->muonPdgId[iL]                = mu.charge*-13;
    looseLeps.push_back(&mu);
    matchLeps.push_back(&mu);
    TVector2 vMu; vMu.SetMagPhi(pt,mu.phi());
    vMETNoMu += vMu;
    gt->nLooseMuon++;
    if (gt->nLooseMuon>=2) break;
  }
  gt->pfmetnomu = vMETNoMu.Mod();

  // now consider all leptons
  gt->nLooseLep = looseLeps.size();
  if (gt->nLooseLep>0) {
    auto ptsort([](panda::Lepton const* l1, panda::Lepton const* l2)->bool {
      return l1->pt() > l2->pt();
    });
    int nToSort = TMath::Min(12,gt->nLooseLep);
    std::partial_sort(looseLeps.begin(),looseLeps.begin()+nToSort,looseLeps.end(),ptsort);
  }
  gt->nTightLep = gt->nTightElectron + gt->nTightMuon;
  if (gt->nLooseLep>0) {
    panda::Lepton* lep1 = looseLeps[0];
    gt->mT = MT(lep1->pt(),lep1->phi(),gt->pfmet,gt->pfmetphi);
  }
  panda::Muon *mu1=0, *mu2=0; panda::Electron *ele1=0, *ele2=0;
  if (gt->nLooseLep>0) { mu1 = dynamic_cast<panda::Muon*>(looseLeps[0]); ele1 = dynamic_cast<panda::Electron*>(looseLeps[0]); };
  if (gt->nLooseLep>1) { mu2 = dynamic_cast<panda::Muon*>(looseLeps[1]); ele2 = dynamic_cast<panda::Electron*>(looseLeps[1]); };
  if (mu1) looseLep1PdgId = mu1->charge*-13; else if(ele1) looseLep1PdgId = ele1->charge*-11;
  if (mu2) looseLep2PdgId = mu2->charge*-13; else if(ele2) looseLep2PdgId = ele2->charge*-11;
  if (gt->nLooseLep>1 && looseLep1PdgId+looseLep2PdgId==0) {
    TLorentzVector v1,v2;
    panda::Lepton *lep1=looseLeps[0], *lep2=looseLeps[1];
    v1.SetPtEtaPhiM(lep1->pt(),lep1->eta(),lep1->phi(),lep1->m());
    v2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());
    gt->diLepMass = (v1+v2).M();
  } else {
    gt->diLepMass = -1;
  }

  tr->TriggerEvent("leptons");
}

void PandaAnalyzer::ComplicatedLeptons() {
  //electrons
  looseLep1PdgId=-1, looseLep2PdgId=-1;
  for (auto& ele : event.electrons) {
    float pt = ele.smearedPt; float eta = ele.eta(); float aeta = fabs(eta);
    if (pt<10 || aeta>2.5 /* || (aeta>1.4442 && aeta<1.566) */) continue;
    if (!ele.veto) continue;
    if (!ElectronIP(ele.eta(),ele.dxy,ele.dz)) continue;
    ele.setPtEtaPhiM(pt,eta,ele.phi(),511e-6);
    unsigned iL=gt->nLooseElectron;
    bool isFake   = ele.hltsafe;
    bool isMedium = ele.medium;
    bool isTight  = ele.tight;
    if (isTight) gt->nTightElectron++;
    int eleSelBit            = kLoose;
    if (isFake  ) eleSelBit |= kFake;
    if (isMedium) eleSelBit |= kMedium;
    if (isTight ) eleSelBit |= kTight;
    gt->electronPt[iL]           = pt;
    gt->electronEta[iL]          = eta;
    gt->electronPhi[iL]          = ele.phi();
    gt->electronD0[iL]           = ele.dxy;
    gt->electronDZ[iL]           = ele.dz;
    gt->electronSfLoose[iL]      = GetCorr(cEleLoose, eta, pt);
    gt->electronSfMedium[iL]     = GetCorr(cEleMedium, eta, pt);
    gt->electronSfTight[iL]      = GetCorr(cEleTight, eta, pt);
    gt->electronSfUnc[iL]        = GetError(cEleMedium, eta, pt);
    gt->electronSfReco[iL]       = GetCorr(cEleReco, eta, pt);
    gt->electronSelBit[iL]       = eleSelBit;
    gt->electronPdgId[iL]        = ele.charge*-11;
    gt->electronChIsoPh[iL]      = ele.chIsoPh;
    gt->electronNhIsoPh[iL]      = ele.nhIsoPh;
    gt->electronPhIsoPh[iL]      = ele.phIsoPh;
    gt->electronEcalIso[iL]      = ele.ecalIso;
    gt->electronHcalIso[iL]      = ele.hcalIso;
    gt->electronTrackIso[iL]     = ele.trackIso;
    gt->electronIsoPUOffset[iL]  = ele.isoPUOffset;
    gt->electronSieie[iL]        = ele.sieie;
    gt->electronSipip[iL]        = ele.sipip;
    gt->electronDEtaInSeed[iL]   = ele.dEtaInSeed;
    gt->electronDPhiIn[iL]       = ele.dPhiIn;
    gt->electronEseed[iL]        = ele.eseed;
    gt->electronHOverE[iL]       = ele.hOverE;
    gt->electronEcalE[iL]        = ele.ecalE;
    gt->electronTrackP[iL]       = ele.trackP;
    gt->electronNMissingHits[iL] = ele.nMissingHits;
    gt->electronTripleCharge[iL] = ele.tripleCharge;
    looseLeps.push_back(&ele);
    matchLeps.push_back(&ele);
    matchEles.push_back(&ele);
    gt->nLooseElectron++;
    if (gt->nLooseElectron>=NLEP) break;

  }

  // muons
  for (auto& mu : event.muons) {
    float pt = mu.pt(); float eta = mu.eta(); float aeta = fabs(eta);
    if (pt<5 || aeta>2.4) continue;
    double ptCorrection=1;
    if (isData) { // perform the rochester correction on the actual particle
      ptCorrection=rochesterCorrection.kScaleDT((int)mu.charge, mu.pt(), mu.eta(), mu.phi(), 0, 0);
    } else { // perform the rochester correction to the simulated particle
      // attempt gen-matching to a final state muon
      bool muonIsTruthMatched=false; TLorentzVector genP4; panda::GenParticle genParticle;
      for (unsigned iG = 0; iG != event.genParticles.size() && !muonIsTruthMatched; ++iG) {
        genParticle = event.genParticles[iG];
        if (genParticle.finalState != 1) continue;
        if (genParticle.pdgid != ((int)mu.charge) * -13) continue;
        genP4.SetPtEtaPhiM(genParticle.pt(), genParticle.eta(), genParticle.phi(), 0.106);
        double dR = genP4.DeltaR(mu.p4());
        if (dR < 0.3) muonIsTruthMatched=true;
      } if (muonIsTruthMatched) { // correct using the gen-particle pt
        double random1=rng.Rndm();
        ptCorrection=rochesterCorrection.kScaleFromGenMC((int)mu.charge, mu.pt(), mu.eta(), mu.phi(), mu.trkLayersWithMmt, genParticle.pt(), random1, 0, 0);
      } else { // if gen match not found, correct the other way
        double random1=rng.Rndm(); double random2=rng.Rndm();
        ptCorrection=rochesterCorrection.kScaleAndSmearMC((int)mu.charge, mu.pt(), mu.eta(), mu.phi(), mu.trkLayersWithMmt, random1, random2, 0, 0);
      }
      pt *= ptCorrection;
    } 
    if (pt<10 || aeta>2.4) continue;
    if (!mu.loose) continue;
    mu.setPtEtaPhiM(pt,eta,mu.phi(),0.106);
    bool isFake   = mu.tight  && mu.combIso()/mu.pt() < 0.4 && mu.chIso/mu.pt() < 0.4;
    bool isMedium = mu.medium && mu.combIso()/mu.pt() < 0.15;
    bool isTight  = mu.tight  && mu.combIso()/mu.pt() < 0.15;
    if (isTight) gt->nTightMuon++;
    int muSelBit            = kLoose;
    if (isFake)   muSelBit |= kFake;
    if (isMedium) muSelBit |= kMedium;
    if (isTight)  muSelBit |= kTight;
    unsigned iL=gt->nLooseMuon;
    gt->muonPt[iL]                   = pt;
    gt->muonEta[iL]                  = eta;
    gt->muonPhi[iL]                  = mu.phi();
    gt->muonD0[iL]                   = mu.dxy;
    gt->muonDZ[iL]                   = mu.dz;
    gt->muonSfLoose[iL]              = GetCorr(cMuLooseID, TMath::Abs(mu.eta()), mu.pt()) * GetCorr(cMuLooseIso, TMath::Abs(mu.eta()), mu.pt());
    gt->muonSfMedium[iL]             = GetCorr(cMuMediumID, TMath::Abs(mu.eta()), mu.pt());
    gt->muonSfTight[iL]              = GetCorr(cMuTightID, TMath::Abs(mu.eta()), mu.pt())  * GetCorr(cMuTightIso, TMath::Abs(mu.eta()), mu.pt());
    gt->muonSfUnc[iL]                = GetError(cMuMediumID , TMath::Abs(mu.eta()), mu.pt());
    gt->muonSfReco[iL]               = GetCorr(cMuReco, mu.eta());
    gt->muonSelBit[iL]               = muSelBit;
    gt->muonPdgId[iL]                = mu.charge*-13;
    gt->muonIsSoftMuon[iL]           = mu.soft;
    gt->muonIsGlobalMuon[iL]         = mu.global;
    gt->muonIsTrackerMuon[iL]        = mu.tracker;
    gt->muonNValidMuon[iL]           = mu.nValidMuon;
    gt->muonNValidPixel[iL]          = mu.nValidPixel;
    gt->muonTrkLayersWithMmt[iL]     = mu.trkLayersWithMmt;
    gt->muonPixLayersWithMmt[iL]     = mu.pixLayersWithMmt;
    gt->muonNMatched[iL]             = mu.nMatched;
    gt->muonChi2LocalPosition[iL]    = mu.chi2LocalPosition;
    gt->muonTrkKink[iL]              = mu.trkKink;
    gt->muonValidFraction[iL]        = mu.validFraction;
    gt->muonNormChi2[iL]             = mu.normChi2;
    gt->muonSegmentCompatibility[iL] = mu.segmentCompatibility;
    looseLeps.push_back(&mu);
    matchLeps.push_back(&mu);
    TVector2 vMu; vMu.SetMagPhi(pt,mu.phi());
    vMETNoMu += vMu;
    gt->nLooseMuon++;
    if (gt->nLooseMuon>=NLEP) break;
  }
  gt->pfmetnomu = vMETNoMu.Mod();

  // now consider all leptons
  gt->nLooseLep = looseLeps.size();
  if (gt->nLooseLep>0) {
    auto ptsort([](panda::Lepton const* l1, panda::Lepton const* l2)->bool {
      return l1->pt() > l2->pt();
    });
    int nToSort = TMath::Min(12,gt->nLooseLep);
    std::partial_sort(looseLeps.begin(),looseLeps.begin()+nToSort,looseLeps.end(),ptsort);
  }
  gt->nTightLep = gt->nTightElectron + gt->nTightMuon;
  if (gt->nLooseLep>0) {
    panda::Lepton* lep1 = looseLeps[0];
    gt->mT = MT(lep1->pt(),lep1->phi(),gt->pfmet,gt->pfmetphi);
  }
  panda::Muon *mu1=0, *mu2=0; panda::Electron *ele1=0, *ele2=0;
  if (gt->nLooseLep>0) { mu1 = dynamic_cast<panda::Muon*>(looseLeps[0]); ele1 = dynamic_cast<panda::Electron*>(looseLeps[0]); };
  if (gt->nLooseLep>1) { mu2 = dynamic_cast<panda::Muon*>(looseLeps[1]); ele2 = dynamic_cast<panda::Electron*>(looseLeps[1]); };
  if (mu1) looseLep1PdgId = mu1->charge*-13; else if(ele1) looseLep1PdgId = ele1->charge*-11;
  if (mu2) looseLep2PdgId = mu2->charge*-13; else if(ele2) looseLep2PdgId = ele2->charge*-11;
  if (gt->nLooseLep>1 && looseLep1PdgId+looseLep2PdgId==0) {
    TLorentzVector v1,v2;
    panda::Lepton *lep1=looseLeps[0], *lep2=looseLeps[1];
    v1.SetPtEtaPhiM(lep1->pt(),lep1->eta(),lep1->phi(),lep1->m());
    v2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());
    gt->diLepMass = (v1+v2).M();
  } else {
    gt->diLepMass = -1;
  }

  tr->TriggerEvent("leptons");
}

void PandaAnalyzer::Photons()
{
    for (auto& pho : event.photons) {
      if (!pho.loose || !pho.csafeVeto)
        continue;
      float pt = pho.pt() * EGMSCALE;
      if (pt<1) continue;
      float eta = pho.eta(), phi = pho.phi();
      if (pt<15 || fabs(eta)>2.5)
        continue;
      loosePhos.push_back(&pho);
      gt->nLoosePhoton++;
      if (gt->nLoosePhoton==1) {
        gt->loosePho1Pt = pt;
        gt->loosePho1Eta = eta;
        gt->loosePho1Phi = phi;
      }
      if ( pho.medium &&
           pt>175 ) { // apply eta cut offline
        if (gt->nLoosePhoton==1)
          gt->loosePho1IsTight=1;
        gt->nTightPhoton++;
        matchPhos.push_back(&pho);
      }
    }

    // TODO - store in a THCorr
    if (isData && gt->nLoosePhoton>0) {
      if (gt->loosePho1Pt>=175 && gt->loosePho1Pt<200)
        gt->sf_phoPurity = 0.04802;
      else if (gt->loosePho1Pt>=200 && gt->loosePho1Pt<250)
        gt->sf_phoPurity = 0.04241;
      else if (gt->loosePho1Pt>=250 && gt->loosePho1Pt<300)
        gt->sf_phoPurity = 0.03641;
      else if (gt->loosePho1Pt>=300 && gt->loosePho1Pt<350)
        gt->sf_phoPurity = 0.0333;
      else if (gt->loosePho1Pt>=350)
        gt->sf_phoPurity = 0.02544;
    }

    tr->TriggerEvent("photons");
}

void PandaAnalyzer::Taus()
{

    for (auto& tau : event.taus) {
      if (analysis->vbf) {
        if (!tau.decayMode || !tau.decayModeNew)
          continue;
        if (!tau.looseIsoMVAOld)
          continue;
      } else {
        if (!tau.decayMode || !tau.decayModeNew)
          continue;
        if (!tau.looseIsoMVA)
          continue;
      }
      if (tau.pt()<18 || fabs(tau.eta())>2.3)
        continue;
      if (IsMatched(&matchLeps,0.16,tau.eta(),tau.phi()))
        continue;
      gt->nTau++;
    }

    tr->TriggerEvent("taus");
}


void PandaAnalyzer::SaveGenLeptons()
{
    gt->genTauPt = -1;
    gt->genElectronPt = -1;
    gt->genMuonPt = -1;
    panda::GenParticle *tau = NULL;
    bool foundTauLeptonic = false; 
    for (auto& gen : event.genParticles) {
      unsigned apdgid = abs(gen.pdgid);
      float pt = gen.pt();
      bool isEmu = false; 

      if (apdgid == 11 && pt > gt->genElectronPt) {
        gt->genElectronPt = pt; 
        gt->genElectronEta = gen.eta(); 
        isEmu = true; 
      }
      
      if (apdgid == 13 && pt > gt->genMuonPt) {
        gt->genMuonPt = pt; 
        gt->genMuonEta = gen.eta(); 
        isEmu = true; 
      }

      if (isEmu && !foundTauLeptonic && tau) {
        const panda::GenParticle *parent = &gen;
        while (parent->parent.isValid()) {
          parent = parent->parent.get();
          if (parent == tau) {
            foundTauLeptonic = true; 
            gt->genTauPt = -1; 
            gt->genTauEta = -1;
            break;
          }
        }
      }

      if (!foundTauLeptonic && apdgid == 15 && pt > gt->genTauPt
          && ((gen.statusFlags & (1 << panda::GenParticle::kIsHardProcess)) != 0 
              || (gen.statusFlags & (1 << panda::GenParticle::kFromHardProcessBeforeFSR)) != 0 
              || ((gen.statusFlags & (1 << panda::GenParticle::kIsDecayedLeptonHadron)) != 0 
                  && (gen.statusFlags & (1 << panda::GenParticle::kFromHardProcess)) != 0
                  )
              )
          ) 
      {
        gt->genTauPt = pt; 
        gt->genTauEta = gen.eta();
      }
    }

    tr->TriggerEvent("gen leptons");

}

void PandaAnalyzer::PhotonSFs()
{
      if (gt->nLoosePhoton < 1)
        return;
      float pt = gt->loosePho1Pt, eta = gt->loosePho1Eta;
      if (gt->loosePho1IsTight)
        gt->sf_pho = GetCorr(cPho,eta,pt);
    tr->TriggerEvent("photon SFs");
}

