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

    if (gt->nLooseElectron>0 && abs(gt->looseLep1PdgId)==11
        && gt->looseLep1IsTight==1) {
      float eff1=0, eff2=0;
      eff1 = GetCorr(cTrigEle,gt->looseLep1Eta,gt->looseLep1Pt);
      if (gt->nLooseElectron>1 && abs(gt->looseLep2PdgId)==11) {
        eff2 = GetCorr(cTrigEle,gt->looseLep2Eta,gt->looseLep2Pt);
      }
      gt->sf_eleTrig = 1 - (1-eff1)*(1-eff2);
    } // done with ele trig SF

    if (gt->nLoosePhoton>0 && gt->loosePho1IsTight)
      gt->sf_phoTrig = GetCorr(cTrigPho,gt->loosePho1Pt);

    if (analysis->vbf) {
      gt->sf_metTrigVBF = GetCorr(cVBF_TrigMET,gt->barrelHTMiss);
      gt->sf_metTrigZmmVBF = GetCorr(cVBF_TrigMETZmm,gt->barrelHTMiss);
    }
    tr->TriggerEvent("triggers");
}

