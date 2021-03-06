#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define EGMSCALE 1

using namespace panda;
using namespace std;


void PandaAnalyzer::SetupJES()
{
  if (uncReader==0) {
    if (isData) {
      TString thisEra = eras.getEra(gt->runNumber);
      for (auto &iter : ak8UncReader) {
        if (! iter.first.Contains("data"))
          continue;
        if (iter.first.Contains(thisEra)) {
          uncReader = iter.second;
          uncReaderAK4 = ak4UncReader[iter.first];
          scaleReaderAK4 = ak4ScaleReader[iter.first];
          break;
        }
      }
    } else {
      uncReader = ak8UncReader["MC"];
      uncReaderAK4 = ak4UncReader["MC"];
      scaleReaderAK4 = ak4ScaleReader["MC"];
    }
  }
}

void PandaAnalyzer::JetBasics() 
{
  gt->barrelJet12Pt = 0;
  gt->barrelHT = 0;
  unsigned nBarrelJets = 0;
  panda::Jet *jet1=0, *jet2=0;
  jot1=0; jot2=0;
  jotUp1=0; jotUp2=0;
  jotDown1=0; jotDown2=0;
  jetUp1=0; jetUp2=0;
  jetDown1=0; jetDown2=0;
  gt->dphipuppimet=999; gt->dphipfmet=999;
  gt->dphipuppiUW=999; gt->dphipfUW=999;
  gt->dphipuppiUZ=999; gt->dphipfUZ=999;
  gt->dphipuppiUA=999; gt->dphipfUA=999;
  float maxJetEta = (analysis->vbf) ? 4.7 : 4.5;
  unsigned nJetDPhi = (analysis->vbf) ? 4 : 5;

  gt->badECALFilter = 1;
  for (auto& jet : *jets) {

    // only do eta-phi checks here
    if (abs(jet.eta()) > maxJetEta)
      continue;
    // NOTE:
    // For VBF we require nTightLep>0, but in monotop looseLep1IsTight
    // No good reason to do that, should switch to former
    // Should update jet cleaning accordingly (just check all loose objects)
    if (IsMatched(&matchLeps,0.16,jet.eta(),jet.phi()) ||
        IsMatched(&matchPhos,0.16,jet.eta(),jet.phi()))
      continue;
    if (analysis->vbf && !jet.loose)
      continue;

    if (analysis->vbf && jet.pt()>20 && fabs(jet.eta())<2.4 && jet.csv>0.8484) {
      ++(gt->jetNMBtags);
    }

    if (jet.pt()>jetPtThreshold) { // nominal jets
      cleanedJets.push_back(&jet);
      if (cleanedJets.size()<3) {
        bool isBad = GetCorr(cBadECALJets,jet.eta(),jet.phi()) > 0;
        if (isBad)
          gt->badECALFilter = 0;
      }

      float csv = (fabs(jet.eta())<2.5) ? jet.csv : -1;
      float cmva = (fabs(jet.eta())<2.5) ? jet.cmva : -1;
      if (fabs(jet.eta())<2.4) {
        centralJets.push_back(&jet);
        if (centralJets.size()==1) {
          jet1 = &jet;
          gt->jet1Pt = jet.pt();
          gt->jet1Eta = jet.eta();
          gt->jet1Phi = jet.phi();
          gt->jet1CSV = csv;
          gt->jet1CMVA = cmva;
          gt->jet1IsTight = jet.monojet ? 1 : 0;
        } else if (centralJets.size()==2) {
          jet2 = &jet;
          gt->jet2Pt = jet.pt();
          gt->jet2Eta = jet.eta();
          gt->jet2Phi = jet.phi();
          gt->jet2CSV = csv;
          gt->jet2CMVA = cmva;
        }
      }

      vJet.SetPtEtaPhiM(jet.pt(),jet.eta(),jet.phi(),jet.m());

      if (analysis->vbf)
        JetVBFBasics(jet);

      if (analysis->monoh || analysis->hbb) {
        JetHbbBasics(jet);
        if (analysis->bjetRegression)
          JetBRegressionInfo(jet);
      }

      // compute dphi wrt mets
      if (cleanedJets.size() <= nJetDPhi) {
        gt->dphipuppimet = std::min(fabs(vJet.DeltaPhi(vPuppiMET)),(double)gt->dphipuppimet);
        gt->dphipfmet = std::min(fabs(vJet.DeltaPhi(vPFMET)),(double)gt->dphipfmet);
        if (analysis->recoil) {
          gt->dphipuppiUA = std::min(fabs(vJet.DeltaPhi(vpuppiUA)),(double)gt->dphipuppiUA);
          gt->dphipuppiUW = std::min(fabs(vJet.DeltaPhi(vpuppiUW)),(double)gt->dphipuppiUW);
          gt->dphipuppiUZ = std::min(fabs(vJet.DeltaPhi(vpuppiUZ)),(double)gt->dphipuppiUZ);
          gt->dphipfUA = std::min(fabs(vJet.DeltaPhi(vpfUA)),(double)gt->dphipfUA);
          gt->dphipfUW = std::min(fabs(vJet.DeltaPhi(vpfUW)),(double)gt->dphipfUW);
          gt->dphipfUZ = std::min(fabs(vJet.DeltaPhi(vpfUZ)),(double)gt->dphipfUZ);
        }
      }
      // btags
      if (csv>0.5426) {
        ++(gt->jetNBtags);
        if (analysis->monoh || analysis->hbb) {
          btaggedJets.push_back(&jet);
          btagindices.push_back(cleanedJets.size()-1);
        }
        if (!analysis->vbf && csv>0.8484) 
          ++(gt->jetNMBtags);
      }
    }

    if (analysis->varyJES)
      JetVaryJES(jet);

  } // VJet loop
  gt->barrelHTMiss = vBarrelJets.Pt();

  switch (gt->whichRecoil) {
    case 0: // MET
      gt->dphipuppiU = gt->dphipuppimet;
      gt->dphipfU = gt->dphipfmet;
      break;
    case -1: // photon
      gt->dphipuppiU = gt->dphipuppiUA;
      gt->dphipfU = gt->dphipfUA;
      break;
    case 1:
      gt->dphipuppiU = gt->dphipuppiUW;
      gt->dphipfU = gt->dphipfUW;
      break;
    case 2:
      gt->dphipuppiU = gt->dphipuppiUZ;
      gt->dphipfU = gt->dphipfUZ;
      break;
    default: // impossible
      break;
  }

  gt->nJet = centralJets.size();
  gt->nJot = cleanedJets.size();

  if (analysis->vbf) {
    JetVBFSystem();
  }

  tr->TriggerEvent("jets");

}

void PandaAnalyzer::JetHbbBasics(panda::Jet& jet)
{
  float csv = (fabs(jet.eta())<2.5) ? jet.csv : -1;
  float cmva = (fabs(jet.eta())<2.5) ? jet.cmva : -1;
  unsigned N = cleanedJets.size()-1;
  gt->jetPt[N]=jet.pt();
  gt->jetEta[N]=jet.eta();
  gt->jetPhi[N]=jet.phi();
  gt->jetE[N]=jet.e();
  gt->jetCSV[N]=csv;
  gt->jetCMVA[N]=cmva;
  gt->jetQGL[N]=jet.qgl;

  tr->TriggerSubEvent("H->bb jet");
}

void PandaAnalyzer::JetBRegressionInfo(panda::Jet& jet)
{
  unsigned N = cleanedJets.size()-1;
  gt->jetEMFrac[N] = jet.cef + jet.nef;
  gt->jetHadFrac[N] = jet.chf + jet.nhf;
  gt->jetLeadingLepPt[N] = 0;
  gt->jetLeadingTrkPt[N] = 0;
  gt->jetNLep[N] = 0;
  for (const panda::Ref<panda::PFCand> &c_iter : jet.constituents) {
    if (!c_iter.isValid())
      continue;
    auto *pf = c_iter.get();
    if (pf->q() != 0) {
      float pt = pf->pt();
      gt->jetLeadingTrkPt[N] = max(pt, gt->jetLeadingTrkPt[N]);
      unsigned pdgid = abs(pf->pdgId());
      if (pdgid == 11 || pdgid == 13) {
        gt->jetNLep[N]++;
        if (pt > gt->jetLeadingLepPt[N]) {
          gt->jetLeadingLepPt[N] = pt;
          gt->jetLeadingLepPtRel[N] = pf->p4().Perp(jet.p4().Vect());
          gt->jetLeadingLepDeltaR[N] = sqrt(DeltaR2(pf->eta(), pf->phi(), jet.eta(), jet.phi()));
        }
      }
    }
  }

  auto& vert = jet.secondaryVertex;
  if (vert.isValid()) {
    gt->jetvtxPt[N] = vert->pt();
    gt->jetvtxMass[N] = vert->m();
    gt->jetvtx3Dval[N] = vert->vtx3DVal;
    gt->jetvtx3Derr[N] = vert->vtx3DeVal;
    gt->jetvtxNtrk[N] = vert->ntrk;
  }

  tr->TriggerSubEvent("b-jet reg info");
}

void PandaAnalyzer::JetVBFBasics(panda::Jet& jet)
{
  if (cleanedJets.size()==1) {
    jot1 = &jet;
    gt->jot1Pt = jet.pt();
    gt->jot1Eta = jet.eta();
    gt->jot1Phi = jet.phi();
    if (analysis->vbf && fabs(gt->jot1Eta)<2.4) { // if it's a central jet, must pass ID requirements
      gt->jot1VBFID = jet.monojet ? 1 : 0;
    } else { // if leading jet is not central, leave the event be
      gt->jot1VBFID = 1;
    }
  } else if (cleanedJets.size()==2) {
    jot2 = &jet;
    gt->jot2Pt = jet.pt();
    gt->jot2Eta = jet.eta();
    gt->jot2Phi = jet.phi();
  }

  if (fabs(jet.eta())<3.0) {
    gt->barrelHT += jet.pt();
    vBarrelJets += vJet;
    if (gt->barrelJet1Pt <= 0) {
      gt->barrelJet1Pt = jet.pt();
      gt->barrelJet1Eta = jet.eta();
    }
  }
  tr->TriggerSubEvent("VBF jet");

}

void PandaAnalyzer::IsoJet(panda::Jet& jet) 
{
  float maxIsoEta = (analysis->monoh) ? 4.5 : 2.5;
  bool isIsoJet = ( (gt->nFatjet==0) || 
      (fabs(jet.eta())<maxIsoEta 
       && DeltaR2(gt->fj1Eta,gt->fj1Phi,jet.eta(),jet.phi())>FATJETMATCHDR2) ); 

  if (isIsoJet) {
    isoJets.push_back(&jet);
    float csv = (fabs(jet.eta())<2.5) ? jet.csv : -1;
    if (csv>0.5426)
      ++gt->isojetNBtags;
    if (isoJets.size()==1) {
      gt->isojet1Pt = jet.pt();
      gt->isojet1CSV = jet.csv;
    } else if (isoJets.size()==2) {
      gt->isojet2Pt = jet.pt();
      gt->isojet2CSV = jet.csv;
    }
    if (analysis->monoh)
      gt->jetIso[cleanedJets.size()-1]=1;
  } else {
    if (analysis->monoh)
      gt->jetIso[cleanedJets.size()-1]=0;
  }
  tr->TriggerSubEvent("iso jets");
}

void PandaAnalyzer::JetVaryJES(panda::Jet& jet)
{
  // do jes variation OUTSIDE of nominal jet check 
  if (jet.ptCorrUp>jetPtThreshold) { // Vary the pT up
    // forward and central jets:
    if (jet.ptCorrUp > gt->jot1PtUp) {
      if (jotUp1) {
        jotUp2 = jotUp1;
        gt->jot2PtUp = gt->jot1PtUp;
        gt->jot2EtaUp = gt->jot1EtaUp;
      }
      jotUp1 = &jet;
      gt->jot1PtUp = jet.ptCorrUp;
      gt->jot1EtaUp = jet.eta();
    } else if (jet.ptCorrUp > gt->jot2PtUp) {
      jotUp2 = &jet;
      gt->jot2PtUp = jet.ptCorrUp;
      gt->jot2EtaUp = jet.eta();
    }
    // central only jets:
    if (fabs(jet.eta()) < 2.4) {
      if (jet.ptCorrUp > gt->jet1PtUp) {
        if (jetUp1) {
          jetUp2 = jetUp1;
          gt->jet2PtUp = gt->jet1PtUp;
          gt->jet2EtaUp = gt->jet1EtaUp;
        }
        jetUp1 = &jet;
        gt->jet1PtUp = jet.ptCorrUp;
        gt->jet1EtaUp = jet.eta();
      } else if (jet.ptCorrUp > gt->jet2PtUp) {
        jetUp2 = &jet;
        gt->jet2PtUp = jet.ptCorrUp;
        gt->jet2EtaUp = jet.eta();
      }
    }
  }
  if (jet.ptCorrDown>jetPtThreshold) { // Vary the pT down:
    // forward and central jets:
    if (jet.ptCorrDown > gt->jot1PtDown) {
      if (jotDown1) {
        jotDown2 = jotDown1;
        gt->jot2PtDown = gt->jot1PtDown;
        gt->jot2EtaDown = gt->jot1EtaDown;
      }
      jotDown1 = &jet;
      gt->jot1PtDown = jet.ptCorrDown;
      gt->jot1EtaDown = jet.eta();
    } else if (jet.ptCorrDown > gt->jot2PtDown) {
      jotDown2 = &jet;
      gt->jot2PtDown = jet.ptCorrDown;
      gt->jot2EtaDown = jet.eta();
    }
    // central only jets:
    if (fabs(jet.eta()) < 2.4) {
      if (jet.ptCorrDown > gt->jet1PtDown) {
        if (jetDown1) {
          jetDown2 = jetDown1;
          gt->jet2PtDown = gt->jet1PtDown;
          gt->jet2EtaDown = gt->jet1EtaDown;
        }
        jetDown1 = &jet;
        gt->jet1PtDown = jet.ptCorrDown;
        gt->jet1EtaDown = jet.eta();
      } else if (jet.ptCorrDown > gt->jet2PtDown) {
        jetDown2 = &jet;
        gt->jet2PtDown = jet.ptCorrDown;
        gt->jet2EtaDown = jet.eta();
      }
    }
  }
  tr->TriggerSubEvent("vary jet JES");
}

void PandaAnalyzer::JetVBFSystem() 
{
  if (gt->nJot>1 && analysis->vbf) {
    TLorentzVector vj1, vj2;
    vj1.SetPtEtaPhiM(jot1->pt(),jot1->eta(),jot1->phi(),jot1->m());
    vj2.SetPtEtaPhiM(jot2->pt(),jot2->eta(),jot2->phi(),jot2->m());
    gt->jot12Mass = (vj1+vj2).M();
    gt->jot12DPhi = vj1.DeltaPhi(vj2);
    gt->jot12DEta = fabs(jot1->eta()-jot2->eta());

    if (analysis->varyJES && jotUp1 && jotUp2) {
      vj1.SetPtEtaPhiM(jotUp1->ptCorrUp,jotUp1->eta(),jotUp1->phi(),jotUp1->m());
      vj2.SetPtEtaPhiM(jotUp2->ptCorrUp,jotUp2->eta(),jotUp2->phi(),jotUp2->m());
      gt->jot12MassUp = (vj1+vj2).M();
      gt->jot12DPhiUp = vj1.DeltaPhi(vj2);
      gt->jot12DEtaUp = fabs(jotUp1->eta()-jotUp2->eta());
    }

    if (analysis->varyJES && jotDown1 && jotDown2) {
      vj1.SetPtEtaPhiM(jotDown1->ptCorrDown,jotDown1->eta(),jotDown1->phi(),jotDown1->m());
      vj2.SetPtEtaPhiM(jotDown2->ptCorrDown,jotDown2->eta(),jotDown2->phi(),jotDown2->m());
      gt->jot12MassDown = (vj1+vj2).M();
      gt->jot12DPhiDown = vj1.DeltaPhi(vj2);
      gt->jot12DEtaDown = fabs(jotDown1->eta()-jotDown2->eta());
    }
  }

  tr->TriggerSubEvent("VBF jet system");
}

void PandaAnalyzer::JetHbbReco() 
{
  float tmp_hbbpt=-99;
  float tmp_hbbeta=-99;
  float tmp_hbbphi=-99;
  float tmp_hbbm=-99;
  int tmp_hbbjtidx1=-1;
  int tmp_hbbjtidx2=-1;
  if (centralJets.size() > 1) {
    vector<Jet*> btagSortedJets = centralJets;
    sort(
      btagSortedJets.begin(),
      btagSortedJets.end(),
      analysis->useCMVA?
        [](panda::Jet *x, panda::Jet *y) -> bool { return x->cmva > y->cmva; } :
        [](panda::Jet *x, panda::Jet *y) -> bool { return x->csv  > y->csv ; }
    );
    map<Jet*, unsigned> order;
    for (unsigned i = 0; i != cleanedJets.size(); ++i) 
      order[cleanedJets[i]] = i;

    panda::Jet *jet_1 = btagSortedJets.at(0);
    TLorentzVector hbbdaughter1;
    hbbdaughter1.SetPtEtaPhiM(jet_1->pt(),jet_1->eta(),jet_1->phi(),jet_1->m());

    panda::Jet *jet_2 = btagSortedJets.at(1);
    TLorentzVector hbbdaughter2;
    hbbdaughter2.SetPtEtaPhiM(jet_2->pt(),jet_2->eta(),jet_2->phi(),jet_2->m());

    TLorentzVector hbbsystem = hbbdaughter1 + hbbdaughter2;

    tmp_hbbpt = hbbsystem.Pt();
    tmp_hbbeta = hbbsystem.Eta();
    tmp_hbbphi = hbbsystem.Phi();
    tmp_hbbm = hbbsystem.M();
    tmp_hbbjtidx1 = order[jet_1];
    tmp_hbbjtidx2 = order[jet_2];
    
  }
  gt->hbbpt = tmp_hbbpt;
  gt->hbbeta = tmp_hbbeta;
  gt->hbbphi = tmp_hbbphi;
  gt->hbbm = tmp_hbbm;
  gt->hbbjtidx[0] = tmp_hbbjtidx1;
  gt->hbbjtidx[1] = tmp_hbbjtidx2;

  
  if (analysis->bjetRegression && gt->hbbm>0.) {
    TLorentzVector hbbdaughters_corr[2];
    
    for (unsigned i = 0; i<2; i++) {
      bjetreg_vars[0] = gt->jetPt[gt->hbbjtidx[i]];
      bjetreg_vars[1] = gt->nJot;
      bjetreg_vars[2] = gt->jetEta[gt->hbbjtidx[i]];
      bjetreg_vars[3] = gt->jetE[gt->hbbjtidx[i]];
      bjetreg_vars[4] = gt->npv;
      bjetreg_vars[5] = gt->jetLeadingTrkPt[gt->hbbjtidx[i]];
      bjetreg_vars[6] = gt->jetLeadingLepPt[gt->hbbjtidx[i]];
      bjetreg_vars[7] = gt->jetNLep[gt->hbbjtidx[i]];
      bjetreg_vars[8] = gt->jetEMFrac[gt->hbbjtidx[i]];
      bjetreg_vars[9] = gt->jetHadFrac[gt->hbbjtidx[i]];
      
      gt->jetRegFac[i] = (bjetreg_reader->EvaluateRegression("BDT method"))[0];
      hbbdaughters_corr[i].SetPtEtaPhiE(gt->jetRegFac[i]*gt->jetPt[gt->hbbjtidx[i]],gt->jetEta[gt->hbbjtidx[i]],gt->jetPhi[gt->hbbjtidx[i]],gt->jetRegFac[i]*gt->jetE[gt->hbbjtidx[i]]);
    }

    TLorentzVector hbbsystem_corr = hbbdaughters_corr[0] + hbbdaughters_corr[1];
    gt->hbbm_reg = hbbsystem_corr.M();
    gt->hbbpt_reg = hbbsystem_corr.Pt();
  }
  
  tr->TriggerEvent("monohiggs");
}

void PandaAnalyzer::GenJetsNu()
{

  std::vector<fastjet::PseudoJet> finalStates;
  std::vector<panda::GenParticle*> bcs;
  for (auto &p : event.genParticles) {
    if (p.finalState && p.pt() > 0.001) {
      finalStates.emplace_back(p.px(), p.py(), p.pz(), p.e());
      continue;
    }
    unsigned apdgid = abs(p.pdgid);
    if (apdgid == 4 || apdgid == 5) {
      bcs.push_back(&p);
    }
  }

  fastjet::ClusterSequenceArea seq(finalStates, *jetDefGen, *areaDef);
  std::vector<fastjet::PseudoJet> allJets(seq.inclusive_jets(0.01));

  genJetsNu.reserve(allJets.size());
  for (auto &pj : allJets) {
    genJetsNu.emplace_back(panda::GenJet());
    genJetsNu.back().setXYZE(pj.px(), pj.py(), pj.pz(), pj.e());
    int flavor = 0;
    for (auto *bc : bcs) {
      if (DeltaR2(pj.eta(), pj.phi(), bc->eta(), bc->phi()) < 0.09) {
        flavor = abs(bc->pdgid);
        break;
      }
    }
    genJetsNu.back().pdgid = flavor;
  }
  tr->TriggerEvent("gen jets nu");
}

