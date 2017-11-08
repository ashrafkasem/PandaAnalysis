#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>
#include "TRandom3.h"
#include "PandaAnalysis/Utilities/src/RoccoR.cc"

#define EGMSCALE 1

using namespace panda;
using namespace std;


void PandaAnalyzer::SimpleLeptons() 
{
    //electrons
    for (auto& ele : event.electrons) {
     float pt = ele.pt()*EGMSCALE; float eta = ele.eta(); float aeta = fabs(eta);
      if (pt<10 || aeta>2.5 /* || (aeta>1.4442 && aeta<1.566) */)
        continue;
      if (!ele.veto)
        continue;
      if (!ElectronIP(ele.eta(),ele.dxy,ele.dz))
        continue;
      looseLeps.push_back(&ele);
      if (analysis->vbf) {
        matchLeps.push_back(&ele);
        matchEles.push_back(&ele);
      }
      gt->nLooseElectron++;
    }

    // muons
    for (auto& mu : event.muons) {
     float pt = mu.pt(); float eta = mu.eta(); float aeta = fabs(eta);
      if (pt<10 || aeta>2.4)
        continue;
      if (!mu.loose)
        continue;
      if (!MuonIsolation(pt,eta,mu.combIso(),panda::kLoose))
        continue;
      looseLeps.push_back(&mu);
      if (analysis->vbf)
        matchLeps.push_back(&mu);
      gt->nLooseMuon++;
      TVector2 vMu; vMu.SetMagPhi(pt,mu.phi());
      vMETNoMu += vMu;
    }
    gt->pfmetnomu = vMETNoMu.Mod();

    // now consider all leptons
    gt->nLooseLep = looseLeps.size();
    if (gt->nLooseLep>0) {
     auto ptsort([](panda::Lepton const* l1, panda::Lepton const* l2)->bool {
       return l1->pt() > l2->pt();
      });
     int nToSort = TMath::Min(3,gt->nLooseLep);
     std::partial_sort(looseLeps.begin(),looseLeps.begin()+nToSort,looseLeps.end(),ptsort);
    }
    int lep_counter=1;
    for (auto* lep : looseLeps) {
      if (lep_counter==1) {
       gt->looseLep1Pt = lep->pt();
       gt->looseLep1Eta = lep->eta();
       gt->looseLep1Phi = lep->phi();
      } else if (lep_counter==2) {
       gt->looseLep2Pt = lep->pt();
       gt->looseLep2Eta = lep->eta();
       gt->looseLep2Phi = lep->phi();
      } else {
        break;
      }
      // now specialize lepton types
      panda::Muon *mu = dynamic_cast<panda::Muon*>(lep);
      if (mu!=NULL) {
        bool isTight = ( mu->tight &&
                MuonIsolation(mu->pt(),mu->eta(),mu->combIso(),panda::kTight) &&
                mu->pt()>20 && fabs(mu->eta())<2.4 );
        if (lep_counter==1) {
          gt->looseLep1PdgId = mu->charge*-13;
          gt->looseLep1IsHLTSafe = 1;
          if (isTight) {
            gt->nTightMuon++;
            gt->looseLep1IsTight = 1;
            if (!analysis->vbf)
              matchLeps.push_back(lep);
          }
        } else if (lep_counter==2) {
          gt->looseLep2PdgId = mu->charge*-13;
          gt->looseLep2IsHLTSafe = 1;
          if (isTight) {
            gt->nTightMuon++;
            gt->looseLep2IsTight = 1;
          }
          if (!analysis->vbf && (isTight || gt->looseLep1IsTight))
            matchLeps.push_back(lep);
        }
      } else {
        panda::Electron *ele = dynamic_cast<panda::Electron*>(lep);
        bool isTight = ( ele->tight &&
                ele->pt()>40 && fabs(ele->eta())<2.5 );
        if (lep_counter==1) {
          gt->looseLep1Pt *= EGMSCALE;
          gt->looseLep1PdgId = ele->charge*-11;
          gt->looseLep1IsHLTSafe = ele->hltsafe ? 1 : 0;
          if (isTight) {
            gt->nTightElectron++;
            gt->looseLep1IsTight = 1;
            if (!analysis->vbf) {
              matchLeps.push_back(lep);
              matchEles.push_back(lep);
            }
          }
        } else if (lep_counter==2) {
          gt->looseLep2Pt *= EGMSCALE;
          gt->looseLep2PdgId = ele->charge*-11;
          gt->looseLep2IsHLTSafe = ele->hltsafe ? 1 : 0;
          if (isTight) {
            gt->nTightElectron++;
            gt->looseLep2IsTight = 1;
          }
          if (!analysis->vbf && (isTight || gt->looseLep1IsTight)) {
            matchLeps.push_back(lep);
            matchEles.push_back(lep);
          }
        }
      }
      ++lep_counter;
    }
    gt->nTightLep = gt->nTightElectron + gt->nTightMuon;
    if (gt->nLooseLep>0) {
      panda::Lepton* lep1 = looseLeps[0];
      gt->mT = MT(lep1->pt(),lep1->phi(),gt->pfmet,gt->pfmetphi);
    }
    if (gt->nLooseLep>1 && gt->looseLep1PdgId+gt->looseLep2PdgId==0) {
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
void PandaAnalyzer::ComplicatedLeptons() 
{
    // TO DO: Hard coded to 2016 rochester corrections for now, need to do this in a better way later
    // Initialize the random seed based on Dylan's age in seconds
    TRandom3 rng; {
     std::time_t t = std::time(0);
     unsigned long int time_now = static_cast<unsigned long int>(time(NULL));
     rng.SetSeed(time_now-731178000);
    } RoccoR rochesterCorrection("PandaAnalysis/data/rcdata.2016.v3");

    //electrons
    for (auto& ele : event.electrons) {
     float pt = ele.smearedPt; float eta = ele.eta(); float aeta = fabs(eta);
      if (pt<10 || aeta>2.5 /* || (aeta>1.4442 && aeta<1.566) */) continue;
      if (!ele.veto) continue;
      if (!ElectronIP(ele.eta(),ele.dxy,ele.dz)) continue;
      ele.setPtEtaPhiM(pt,eta,ele.phi(),511e-6);
      looseLeps.push_back(&ele);
      if (analysis->vbf) {
        matchLeps.push_back(&ele);
        matchEles.push_back(&ele);
      }
      gt->nLooseElectron++;
    }

    // muons
    for (auto& mu : event.muons) {
     float pt = mu.pt(); float eta = mu.eta(); float aeta = fabs(eta);
     double ptCorrection=1;
     if(isData) { // perform the rochester correction on the actual particle
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
      } if(muonIsTruthMatched) { // correct using the gen-particle pt
       double random1=rng.Rndm();
       ptCorrection=rochesterCorrection.kScaleFromGenMC((int)mu.charge, mu.pt(), mu.eta(), mu.phi(), mu.trkLayersWithMmt, genParticle.pt(), random1, 0, 0);
      } else { // if gen match not found, correct the other way
       double random1=rng.Rndm(); double random2=rng.Rndm();
       ptCorrection=rochesterCorrection.kScaleAndSmearMC((int)mu.charge, mu.pt(), mu.eta(), mu.phi(), mu.trkLayersWithMmt, random1, random2, 0, 0);
      }
      pt *= ptCorrection;
     } 
     if (pt<10 || aeta>2.4) continue;
     if(!mu.loose) continue;
     mu.setPtEtaPhiM(pt,eta,mu.phi(),0.106);
     looseLeps.push_back(&mu);
     matchLeps.push_back(&mu);
     gt->nLooseMuon++;
     TVector2 vMu; vMu.SetMagPhi(pt,mu.phi());
     vMETNoMu += vMu;
    }
    gt->pfmetnomu = vMETNoMu.Mod();

    // now consider all leptons
    gt->nLooseLep = looseLeps.size();
    if (gt->nLooseLep>0) {
     auto ptsort([](panda::Lepton const* l1, panda::Lepton const* l2)->bool {
       return l1->pt() > l2->pt();
      });
     int nToSort = TMath::Min(3,gt->nLooseLep);
     std::partial_sort(looseLeps.begin(),looseLeps.begin()+nToSort,looseLeps.end(),ptsort);
    }
    unsigned lep_counter=1;
    for (auto* lep : looseLeps) {
      // assign general lepton properties
      if(lep_counter>4) break; 
      if      (lep_counter==1){ gt->looseLep1Pt=lep->pt(); gt->looseLep1Eta=lep->eta(); gt->looseLep1Phi=lep->phi();}
      else if (lep_counter==2){ gt->looseLep2Pt=lep->pt(); gt->looseLep2Eta=lep->eta(); gt->looseLep2Phi=lep->phi();}
      else if (lep_counter==3){ gt->looseLep3Pt=lep->pt(); gt->looseLep3Eta=lep->eta(); gt->looseLep3Phi=lep->phi();}
      else if (lep_counter==4){ gt->looseLep4Pt=lep->pt(); gt->looseLep4Eta=lep->eta(); gt->looseLep4Phi=lep->phi();}
      else break;
      // now specialize lepton types
      panda::Muon *mu = dynamic_cast<panda::Muon*>(lep); if (mu!=NULL) {
       // assign muon specific properties
       int *lepPdgId, *isHLTSafe, *lepSelBit; 
       float *lepPt, *lepEta, *lepPhi, *lepSfTrk, *lepSfLoose, *lepSfMedium, *lepSfTight, *lepSfUnc;
       int *soft, *global, *tracker;
       int *nValidMuon, *nValidPixel, *trkLayersWithMmt, *pixLayersWithMmt, *nMatched, *chi2LocalPosition, *trkKink;
       float *validFraction, *normChi2, *segmentCompatibility;
       if (lep_counter==1) {
        lepPdgId             = &(gt->looseLep1PdgId                ); 
        isHLTSafe            = &(gt->looseLep1IsHLTSafe            ); 
        lepSelBit            = &(gt->looseLep1SelBit               ); 
        lepSfTrk             = &(gt->sf_trk1                       ); 
        lepSfLoose           = &(gt->sf_loose1                     ); 
        lepSfMedium          = &(gt->sf_medium1                    ); 
        lepSfTight           = &(gt->sf_tight1                     ); 
        lepSfUnc             = &(gt->sf_unc1                       ); 
        soft                 = &(gt->looseLep1IsSoftMuon           ); 
        global               = &(gt->looseLep1IsGlobalMuon         ); 
        tracker              = &(gt->looseLep1IsTrackerMuon        ); 
        nValidMuon           = &(gt->looseLep1NValidMuon           ); 
        nValidPixel          = &(gt->looseLep1NValidPixel          ); 
        trkLayersWithMmt     = &(gt->looseLep1TrkLayersWithMmt     ); 
        pixLayersWithMmt     = &(gt->looseLep1PixLayersWithMmt     ); 
        nMatched             = &(gt->looseLep1NMatched             ); 
        chi2LocalPosition    = &(gt->looseLep1Chi2LocalPosition    ); 
        trkKink              = &(gt->looseLep1TrkKink              ); 
        validFraction        = &(gt->looseLep1ValidFraction        ); 
        normChi2             = &(gt->looseLep1NormChi2             ); 
        segmentCompatibility = &(gt->looseLep1SegmentCompatibility ); 
       } else if(lep_counter==2) { 
        lepPdgId             = &(gt->looseLep2PdgId                ); 
        isHLTSafe            = &(gt->looseLep2IsHLTSafe            ); 
        lepSelBit            = &(gt->looseLep2SelBit               ); 
        lepSfTrk             = &(gt->sf_trk2                       ); 
        lepSfLoose           = &(gt->sf_loose2                     ); 
        lepSfMedium          = &(gt->sf_medium2                    ); 
        lepSfTight           = &(gt->sf_tight2                     ); 
        lepSfUnc             = &(gt->sf_unc2                       ); 
        soft                 = &(gt->looseLep2IsSoftMuon           ); 
        global               = &(gt->looseLep2IsGlobalMuon         ); 
        tracker              = &(gt->looseLep2IsTrackerMuon        ); 
        nValidMuon           = &(gt->looseLep2NValidMuon           ); 
        nValidPixel          = &(gt->looseLep2NValidPixel          ); 
        trkLayersWithMmt     = &(gt->looseLep2TrkLayersWithMmt     ); 
        pixLayersWithMmt     = &(gt->looseLep2PixLayersWithMmt     ); 
        nMatched             = &(gt->looseLep2NMatched             ); 
        chi2LocalPosition    = &(gt->looseLep2Chi2LocalPosition    ); 
        trkKink              = &(gt->looseLep2TrkKink              ); 
        validFraction        = &(gt->looseLep2ValidFraction        ); 
        normChi2             = &(gt->looseLep2NormChi2             ); 
        segmentCompatibility = &(gt->looseLep2SegmentCompatibility ); 
       } else if(lep_counter==3) { 
        lepPdgId             = &(gt->looseLep3PdgId                ); 
        isHLTSafe            = &(gt->looseLep3IsHLTSafe            ); 
        lepSelBit            = &(gt->looseLep3SelBit               ); 
        lepSfTrk             = &(gt->sf_trk3                       ); 
        lepSfLoose           = &(gt->sf_loose3                     ); 
        lepSfMedium          = &(gt->sf_medium3                    ); 
        lepSfTight           = &(gt->sf_tight3                     ); 
        lepSfUnc             = &(gt->sf_unc3                       ); 
        soft                 = &(gt->looseLep3IsSoftMuon           ); 
        global               = &(gt->looseLep3IsGlobalMuon         ); 
        tracker              = &(gt->looseLep3IsTrackerMuon        ); 
        nValidMuon           = &(gt->looseLep3NValidMuon           ); 
        nValidPixel          = &(gt->looseLep3NValidPixel          ); 
        trkLayersWithMmt     = &(gt->looseLep3TrkLayersWithMmt     ); 
        pixLayersWithMmt     = &(gt->looseLep3PixLayersWithMmt     ); 
        nMatched             = &(gt->looseLep3NMatched             ); 
        chi2LocalPosition    = &(gt->looseLep3Chi2LocalPosition    ); 
        trkKink              = &(gt->looseLep3TrkKink              ); 
        validFraction        = &(gt->looseLep3ValidFraction        ); 
        normChi2             = &(gt->looseLep3NormChi2             ); 
        segmentCompatibility = &(gt->looseLep3SegmentCompatibility ); 
       } else if(lep_counter==4) { 
        lepPdgId             = &(gt->looseLep4PdgId                ); 
        isHLTSafe            = &(gt->looseLep4IsHLTSafe            ); 
        lepSelBit            = &(gt->looseLep4SelBit               ); 
        lepSfTrk             = &(gt->sf_trk4                       ); 
        lepSfLoose           = &(gt->sf_loose4                     ); 
        lepSfMedium          = &(gt->sf_medium4                    ); 
        lepSfTight           = &(gt->sf_tight4                     ); 
        lepSfUnc             = &(gt->sf_unc4                       ); 
        soft                 = &(gt->looseLep4IsSoftMuon           ); 
        global               = &(gt->looseLep4IsGlobalMuon         ); 
        tracker              = &(gt->looseLep4IsTrackerMuon        ); 
        nValidMuon           = &(gt->looseLep4NValidMuon           ); 
        nValidPixel          = &(gt->looseLep4NValidPixel          ); 
        trkLayersWithMmt     = &(gt->looseLep4TrkLayersWithMmt     ); 
        pixLayersWithMmt     = &(gt->looseLep4PixLayersWithMmt     ); 
        nMatched             = &(gt->looseLep4NMatched             ); 
        chi2LocalPosition    = &(gt->looseLep4Chi2LocalPosition    ); 
        trkKink              = &(gt->looseLep4TrkKink              ); 
        validFraction        = &(gt->looseLep4ValidFraction        ); 
        normChi2             = &(gt->looseLep4NormChi2             ); 
        segmentCompatibility = &(gt->looseLep4SegmentCompatibility ); 
       }
       *lepPdgId = mu->charge*-13;
       *isHLTSafe=1;
       bool isFake   = mu->tight  && mu->combIso()/mu->pt() < 0.4 && mu->chIso/mu->pt() < 0.4;
       bool isMedium = mu->medium && mu->combIso()/mu->pt() < 0.15;
       bool isTight  = mu->tight  && mu->combIso()/mu->pt() < 0.15;
       *lepSelBit              = kLoose; 
       if(isFake)   *lepSelBit |= kFake; 
       if(isMedium) *lepSelBit |= kMedium; 
       if(isTight)  *lepSelBit |= kTight;
       *lepSfTrk = GetCorr(cMuReco,mu->eta());
       *lepSfLoose = GetCorr(cMuLooseID,TMath::Abs(mu->eta()),mu->pt()) * GetCorr(cMuLooseIso,TMath::Abs(mu->eta()),mu->pt());
       *lepSfMedium = GetCorr(cMuMediumID,TMath::Abs(mu->eta()),mu->pt());
       *lepSfTight = GetCorr(cMuTightID,TMath::Abs(mu->eta()),mu->pt())  * GetCorr(cMuTightIso,TMath::Abs(mu->eta()),mu->pt());
	   *lepSfUnc = GetError(cMuMediumID ,TMath::Abs(mu->eta()),mu->pt());
       *soft                 = mu->soft                ; 
       *global               = mu->global              ; 
       *tracker              = mu->tracker             ; 
       *nValidMuon           = mu->nValidMuon          ; 
       *nValidPixel          = mu->nValidPixel         ; 
       *trkLayersWithMmt     = mu->trkLayersWithMmt    ; 
       *pixLayersWithMmt     = mu->pixLayersWithMmt    ; 
       *nMatched             = mu->nMatched            ; 
       *chi2LocalPosition    = mu->chi2LocalPosition   ; 
       *trkKink              = mu->trkKink             ; 
       *validFraction        = mu->validFraction       ; 
       *normChi2             = mu->normChi2            ; 
       *segmentCompatibility = mu->segmentCompatibility; 
       gt->nTightMuon+=isTight;
      } else {
        panda::Electron *ele = dynamic_cast<panda::Electron*>(lep); assert(ele!=NULL);
        bool isTight = (ele->tight && ele->pt()>10 && fabs(ele->eta())<2.5);
        if (lep_counter==1) {
          gt->looseLep1Pt *= EGMSCALE;
          gt->looseLep1PdgId = ele->charge*-11;
          gt->looseLep1IsHLTSafe = ele->hltsafe ? 1 : 0;
          if (isTight) {
            gt->nTightElectron++;
            gt->looseLep1IsTight = 1;
            if (!analysis->vbf) {
              matchLeps.push_back(lep);
              matchEles.push_back(lep);
            }
          }
        } else if (lep_counter==2) {
          gt->looseLep2Pt *= EGMSCALE;
          gt->looseLep2PdgId = ele->charge*-11;
          gt->looseLep2IsHLTSafe = ele->hltsafe ? 1 : 0;
          if (isTight) {
            gt->nTightElectron++;
            gt->looseLep2IsTight = 1;
          }
          if (!analysis->vbf && (isTight || gt->looseLep1IsTight)) {
            matchLeps.push_back(lep);
            matchEles.push_back(lep);
          }
        }
      }
      ++lep_counter;
    }
    gt->nTightLep = gt->nTightElectron + gt->nTightMuon;
    if (gt->nLooseLep>0) {
      panda::Lepton* lep1 = looseLeps[0];
      gt->mT = MT(lep1->pt(),lep1->phi(),gt->pfmet,gt->pfmetphi);
    }
    if (gt->nLooseLep>1 && gt->looseLep1PdgId+gt->looseLep2PdgId==0) {
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

void PandaAnalyzer::LeptonSFs()
{
      for (unsigned int iL=0; iL!=TMath::Min(gt->nLooseLep,2); ++iL) {
        auto* lep = looseLeps.at(iL);
        float pt = lep->pt(), eta = lep->eta(), aeta = TMath::Abs(eta);
        bool isTight = (iL==0 && gt->looseLep1IsTight) || (iL==1 && gt->looseLep2IsTight);
        auto* mu = dynamic_cast<panda::Muon*>(lep);
        if (mu!=NULL) {
          if (isTight) {
            gt->sf_lepID *= GetCorr(cMuTightID,aeta,pt);
            gt->sf_lepIso *= GetCorr(cMuTightIso,aeta,pt);
          } else {
            gt->sf_lepID *= GetCorr(cMuLooseID,aeta,pt);
            gt->sf_lepIso *= GetCorr(cMuLooseIso,aeta,pt);
          }
          gt->sf_lepTrack *= GetCorr(cMuReco,gt->npv);
        } else {
          if (isTight) {
            gt->sf_lepID *= GetCorr(cEleTight,eta,pt);
          } else {
            gt->sf_lepID *= GetCorr(cEleVeto,eta,pt);
          }
          gt->sf_lepTrack *= GetCorr(cEleReco,eta,pt);
        }
      }

    tr->TriggerEvent("lepton SFs");
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

