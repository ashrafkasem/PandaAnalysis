#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define EGMSCALE 1

using namespace panda;
using namespace std;


void PandaAnalyzer::RegisterTriggers() 
{
  for (auto &th : triggerHandlers) {
    unsigned N = th.paths.size();
    for (unsigned i = 0; i != N; i++) {
      unsigned panda_idx = event.registerTrigger(th.paths.at(i));
      th.indices[i] = panda_idx;
      if (DEBUG) PDebug("PandaAnalyzer::RegisterTriggers",
        Form("Got index %d for trigger path %s", panda_idx, th.paths.at(i).Data())
      );
    }
  }
}

bool PandaAnalyzer::RecoilPresel() 
{
    if ( (preselBits&kMonotop) || (preselBits&kMonohiggs) || 
         (preselBits&kMonojet) || (preselBits&kRecoil) ) 
    {
       if (event.recoil.max<175)
         return false;
    } 
    return true;
}

void PandaAnalyzer::TriggerEffs()
{

    // trigger efficiencies
    gt->sf_metTrig = GetCorr(cTrigMET,gt->pfmetnomu);
    gt->sf_metTrigZmm = GetCorr(cTrigMETZmm,gt->pfmetnomu);

    if (gt->nLooseElectron>0) {
      panda::Electron *ele1=0, *ele2=0;
      if (gt->nLooseLep>0) ele1 = dynamic_cast<panda::Electron*>(looseLeps[0]);
      if (gt->nLooseLep>1) ele2 = dynamic_cast<panda::Electron*>(looseLeps[1]);
      float eff1=0, eff2=0;
      if (ele1 && ele1->tight) {
        eff1 = GetCorr(cTrigEle, ele1->eta(), ele1->pt());
        if (ele2 && ele2->tight)
          eff2 = GetCorr(cTrigEle, ele2->eta(), ele2->pt());
        gt->sf_eleTrig = 1 - (1-eff1)*(1-eff2);
      }
    } // done with ele trig SF

    if (gt->nLoosePhoton>0 && gt->loosePho1IsTight)
      gt->sf_phoTrig = GetCorr(cTrigPho,gt->loosePho1Pt);

    if (analysis->vbf) {
      gt->sf_metTrigVBF = GetCorr(cVBF_TrigMET,gt->barrelHTMiss);
      gt->sf_metTrigZmmVBF = GetCorr(cVBF_TrigMETZmm,gt->barrelHTMiss);
    }
    tr->TriggerEvent("triggers");
}

void PandaAnalyzer::Recoil()
{
    TLorentzVector vpfUp; vpfUp.SetPtEtaPhiM(gt->pfmetUp,0,gt->pfmetphi,0);
    TLorentzVector vpfDown; vpfDown.SetPtEtaPhiM(gt->pfmetDown,0,gt->pfmetphi,0);
    TLorentzVector vpfUWUp, vpfUWDown;
    TLorentzVector vObj1, vObj2;
    gt->whichRecoil = 0; // -1=photon, 0=MET, 1,2=nLep
    if (gt->nLooseLep>0) {
      panda::Lepton *lep1 = looseLeps.at(0);
      vObj1.SetPtEtaPhiM(lep1->pt(),lep1->eta(),lep1->phi(),lep1->m());

      // one lep => W
      vpuppiUW = vPuppiMET+vObj1; gt->puppiUWmag=vpuppiUW.Pt(); gt->puppiUWphi=vpuppiUW.Phi();
      vpfUW = vPFMET+vObj1; gt->pfUWmag=vpfUW.Pt(); gt->pfUWphi=vpfUW.Phi();
      
      if (analysis->varyJES) {
        vpfUWUp = vpfUp+vObj1; gt->pfUWmagUp = vpfUWUp.Pt();
        vpfUWDown = vpfDown+vObj1; gt->pfUWmagDown = vpfUWDown.Pt();
      }

      if (gt->nLooseLep>1 && looseLep1PdgId+looseLep2PdgId==0) {
        // two OS lep => Z
        panda::Lepton *lep2 = looseLeps.at(1);
        vObj2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());

        vpuppiUZ=vpuppiUW+vObj2; gt->puppiUZmag=vpuppiUZ.Pt(); gt->puppiUZphi=vpuppiUZ.Phi();
        vpfUZ=vpfUW+vObj2; gt->pfUZmag=vpfUZ.Pt(); gt->pfUZphi=vpfUZ.Phi();

        if (analysis->varyJES) {
          TLorentzVector vpfUZUp = vpfUWUp+vObj2; gt->pfUZmagUp = vpfUZUp.Pt();
          TLorentzVector vpfUZDown = vpfUWDown+vObj2; gt->pfUZmagDown = vpfUZDown.Pt();
        }

        vpuppiU = vpuppiUZ; vpfU = vpfUZ;
        gt->whichRecoil = 2;
      } else {
        vpuppiU = vpuppiUW; vpfU = vpfUW;
        gt->whichRecoil = 1;
      }
    }
    if (gt->nLoosePhoton>0) {
      panda::Photon *pho = loosePhos.at(0);
      vObj1.SetPtEtaPhiM(pho->pt(),pho->eta(),pho->phi(),0.);

      vpuppiUA=vPuppiMET+vObj1; gt->puppiUAmag=vpuppiUA.Pt(); gt->puppiUAphi=vpuppiUA.Phi();
      vpfUA=vPFMET+vObj1; gt->pfUAmag=vpfUA.Pt(); gt->pfUAphi=vpfUA.Phi();

      if (analysis->varyJES) {
        TLorentzVector vpfUAUp = vpfUp+vObj1; gt->pfUAmagUp = vpfUAUp.Pt();
        TLorentzVector vpfUADown = vpfDown+vObj1; gt->pfUAmagDown = vpfUADown.Pt();
      }

      if (gt->nLooseLep==0) {
        vpuppiU = vpuppiUA; vpfU = vpfUA;
        gt->whichRecoil = -1;
      }
    }
    if (gt->nLooseLep==0 && gt->nLoosePhoton==0) {
      vpuppiU = vPuppiMET;
      vpfU = vPFMET;
      gt->whichRecoil = 0;
    }
    gt->puppiUmag = vpuppiU.Pt();
    gt->puppiUphi = vpuppiU.Phi();
    gt->pfUmag = vpfU.Pt();
    gt->pfUphi = vpfU.Phi();

    tr->TriggerEvent("recoils");
}

void PandaAnalyzer::GenPartonStudy() 
{
  // For now, simple B and C counting
  for (auto& gen : event.genParticles) {
    float pt = gen.pt();
    int pdgid = gen.pdgid;
    if (gen.parent.isValid() && gen.parent->pdgid==gen.pdgid)
      continue;
    //count bs and cs
    int apdgid = abs(pdgid);
    if (apdgid!=5 && apdgid!=4) 
      continue;
    if (gen.pt()>5) {
      gt->nHF++;
      if (apdgid==5)
        gt->nB++;
    }
  }
}

void PandaAnalyzer::GetMETSignificance()
{
  float pfEt = 0;
  float puppiEt = 0;

  for (auto& pfCand : event.pfCandidates){
    TLorentzVector pfcand(0,0,0,0);
    pfcand.SetPtEtaPhiM(pfCand.pt(),pfCand.eta(),pfCand.phi(),pfCand.m());
    puppiEt += pfcand.Et()*pfCand.puppiW();
    pfEt += pfcand.Et();
  }

  gt->pfmetsig = event.pfMet.pt/sqrt(pfEt);
  gt->puppimetsig = event.puppiMet.pt/sqrt(puppiEt);
}

