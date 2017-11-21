#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define EGMSCALE 1

using namespace panda;
using namespace std;

void PandaAnalyzer::CalcBJetSFs(BTagType bt, int flavor,
                                double eta, double pt, double eff, double uncFactor,
                                double &sf, double &sfUp, double &sfDown) 
{
  if (flavor==5) {
    sf     = btagReaders[bt]->eval_auto_bounds("central",BTagEntry::FLAV_B,eta,pt);
    sfUp   = btagReaders[bt]->eval_auto_bounds("up",BTagEntry::FLAV_B,eta,pt);
    sfDown = btagReaders[bt]->eval_auto_bounds("down",BTagEntry::FLAV_B,eta,pt);
  } else if (flavor==4) {
    sf     = btagReaders[bt]->eval_auto_bounds("central",BTagEntry::FLAV_C,eta,pt);
    sfUp   = btagReaders[bt]->eval_auto_bounds("up",BTagEntry::FLAV_C,eta,pt);
    sfDown = btagReaders[bt]->eval_auto_bounds("down",BTagEntry::FLAV_C,eta,pt);
  } else {
    sf     = btagReaders[bt]->eval_auto_bounds("central",BTagEntry::FLAV_UDSG,eta,pt);
    sfUp   = btagReaders[bt]->eval_auto_bounds("up",BTagEntry::FLAV_UDSG,eta,pt);
    sfDown = btagReaders[bt]->eval_auto_bounds("down",BTagEntry::FLAV_UDSG,eta,pt);
  }

  sfUp = uncFactor*(sfUp-sf)+sf;
  sfDown = uncFactor*(sfDown-sf)+sf;
  return;
}

void PandaAnalyzer::EvalBTagSF(std::vector<btagcand> &cands, std::vector<double> &sfs,
                               GeneralTree::BTagShift shift,GeneralTree::BTagJet jettype, bool do2) 
{
  float sf0 = 1, sf1 = 1, sfGT0 = 1, sf2=1;
  float prob_mc0=1, prob_data0=1;
  float prob_mc1=0, prob_data1=0;
  unsigned int nC = cands.size();

  for (unsigned int iC=0; iC!=nC; ++iC) {
    double sf_i = sfs[iC];
    double eff_i = cands[iC].eff;
    prob_mc0 *= (1-eff_i);
    prob_data0 *= (1-sf_i*eff_i);
    float tmp_mc1=1, tmp_data1=1;
    for (unsigned int jC=0; jC!=nC; ++jC) {
      if (iC==jC) continue;
      double sf_j = sfs[jC];
      double eff_j = cands[jC].eff;
      tmp_mc1 *= (1-eff_j);
      tmp_data1 *= (1-eff_j*sf_j);
    }
    prob_mc1 += eff_i * tmp_mc1;
    prob_data1 += eff_i * sf_i * tmp_data1;
  }
  
  if (nC>0) {
    sf0 = prob_data0/prob_mc0;
    sf1 = prob_data1/prob_mc1;
    sfGT0 = (1-prob_data0)/(1-prob_mc0);
  }

  GeneralTree::BTagParams p;
  p.shift = shift;
  p.jet = jettype;
  p.tag=GeneralTree::b0; gt->sf_btags[p] = sf0;
  p.tag=GeneralTree::b1; gt->sf_btags[p] = sf1;
  p.tag=GeneralTree::bGT0; gt->sf_btags[p] = sfGT0;

  if (do2) {
    float prob_mc2=0, prob_data2=0;
    unsigned int nC = cands.size();

    for (unsigned int iC=0; iC!=nC; ++iC) {
      double sf_i = sfs[iC], eff_i = cands[iC].eff;
      for (unsigned int jC=iC+1; jC!=nC; ++jC) {
        double sf_j = sfs[jC], eff_j = cands[jC].eff;
        float tmp_mc2=1, tmp_data2=1;
        for (unsigned int kC=0; kC!=nC; ++kC) {
          if (kC==iC || kC==jC) continue;
          double sf_k = sfs[kC], eff_k = cands[kC].eff;
          tmp_mc2 *= (1-eff_k);
          tmp_data2 *= (1-eff_k*sf_k);
        }
        prob_mc2 += eff_i * eff_j * tmp_mc2;
        prob_data2 += eff_i * sf_i * eff_j * sf_j * tmp_data2;
      }
    }

    if (nC>1) {
      sf2 = prob_data2/prob_mc2;
    }

    p.tag=GeneralTree::b2; gt->sf_btags[p] = sf2;
  }

}


void PandaAnalyzer::JetBtagSFs() 
{
      // now get the jet btag SFs
      vector<btagcand> btagcands;
      vector<double> sf_cent, sf_bUp, sf_bDown, sf_mUp, sf_mDown;

      unsigned int nJ = centralJets.size();
      for (unsigned int iJ=0; iJ!=nJ; ++iJ) {
        panda::Jet *jet = centralJets.at(iJ);
        bool isIsoJet=false;
        if (std::find(isoJets.begin(), isoJets.end(), jet) != isoJets.end())
          isIsoJet = true;
        int flavor=0;
        float genpt=0;
        for (auto& gen : event.genParticles) {
          int apdgid = abs(gen.pdgid);
          if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
            continue;
          double dr2 = DeltaR2(jet->eta(),jet->phi(),gen.eta(),gen.phi());
          if (dr2<0.09) {
            genpt = gen.pt();
            if (apdgid==4 || apdgid==5) {
              flavor=apdgid;
              break;
            } else {
              flavor=0;
            }
          }
        } // finding the jet flavor
        float pt = jet->pt();
        float btagUncFactor = 1;
        float eta = jet->eta();
        double eff(1),sf(1),sfUp(1),sfDown(1);
        unsigned int binpt = btagpt.bin(pt);
        unsigned int bineta = btageta.bin(fabs(eta));
        if (flavor==5)
          eff = beff[bineta][binpt];
        else if (flavor==4)
          eff = ceff[bineta][binpt];
        else
          eff = lfeff[bineta][binpt];
        if (jet==centralJets.at(0)) {
          gt->jet1Flav = flavor;
          gt->jet1GenPt = genpt;
        } else if (jet==centralJets.at(1)) {
          gt->jet2Flav = flavor;
          gt->jet2GenPt = genpt;
        }
        if (isIsoJet) {
          if (jet==isoJets.at(0))
            gt->isojet1Flav = flavor;
          else if (jet==isoJets.at(1))
            gt->isojet2Flav = flavor;

          CalcBJetSFs(bJetL,flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
          btagcands.push_back(btagcand(iJ,flavor,eff,sf,sfUp,sfDown));
          sf_cent.push_back(sf);

          if (flavor>0) {
            sf_bUp.push_back(sfUp); sf_bDown.push_back(sfDown);
            sf_mUp.push_back(sf); sf_mDown.push_back(sf);
          } else {
            sf_bUp.push_back(sf); sf_bDown.push_back(sf);
            sf_mUp.push_back(sfUp); sf_mDown.push_back(sfDown);
          }

        }

      } // loop over jets

      EvalBTagSF(btagcands,sf_cent,GeneralTree::bCent,GeneralTree::bJet);
      EvalBTagSF(btagcands,sf_bUp,GeneralTree::bBUp,GeneralTree::bJet);
      EvalBTagSF(btagcands,sf_bDown,GeneralTree::bBDown,GeneralTree::bJet);
      EvalBTagSF(btagcands,sf_mUp,GeneralTree::bMUp,GeneralTree::bJet);
      EvalBTagSF(btagcands,sf_mDown,GeneralTree::bMDown,GeneralTree::bJet);

    tr->TriggerEvent("ak4 gen-matching");
}

void PandaAnalyzer::JetCMVAWeights() 
{
   if (analysis->useCMVA) {
    gt->sf_cmvaweight = 1;
    gt->sf_cmva_LFup = 1;
    gt->sf_cmva_LFdown = 1;
    gt->sf_cmva_HFup = 1;
    gt->sf_cmva_HFdown = 1;
    gt->sf_cmva_HFStats1up = 1;
    gt->sf_cmva_HFStats1down = 1;
    gt->sf_cmva_LFStats1up = 1;
    gt->sf_cmva_LFStats1down = 1;
    gt->sf_cmva_HFStats2up = 1;
    gt->sf_cmva_HFStats2down = 1;
    gt->sf_cmva_LFStats2up = 1;
    gt->sf_cmva_LFStats2down = 1;
    gt->sf_cmva_CErr1up = 1;
    gt->sf_cmva_CErr1down = 1;
    gt->sf_cmva_CErr2up = 1;
    gt->sf_cmva_CErr2down = 1;
  } else {
    gt->sf_csvweight = 1;
    gt->sf_csv_LFup = 1;
    gt->sf_csv_LFdown = 1;
    gt->sf_csv_HFup = 1;
    gt->sf_csv_HFdown = 1;
    gt->sf_csv_HFStats1up = 1;
    gt->sf_csv_HFStats1down = 1;
    gt->sf_csv_LFStats1up = 1;
    gt->sf_csv_LFStats1down = 1;
    gt->sf_csv_HFStats2up = 1;
    gt->sf_csv_HFStats2down = 1;
    gt->sf_csv_LFStats2up = 1;
    gt->sf_csv_LFStats2down = 1;
    gt->sf_csv_CErr1up = 1;
    gt->sf_csv_CErr1down = 1;
    gt->sf_csv_CErr2up = 1;
    gt->sf_csv_CErr2down = 1;
  }
  if (centralJets.size() < 1) return;

  //get vectors of jet properties
  std::vector<double> jetPts, jetEtas, jetCSVs, jetCMVAs;
  std::vector<int> jetFlavors;
  for (unsigned i = 0; i != centralJets.size(); ++i) {
    panda::Jet* jet=centralJets[i];
    jetPts.push_back(jet->pt());
    jetEtas.push_back(jet->eta());
    jetCSVs.push_back(jet->csv); 
    jetCMVAs.push_back(jet->cmva);
    int flavor=0; // Will be considered a LF jet if we can't match to a gen particle
    for (auto &gen : event.ak4GenJets) {
      if (DeltaR2(gen.eta(), gen.phi(), jet->eta(), jet->phi()) < 0.09) {
        flavor=gen.pdgid;
        break;
      }
    }
    jetFlavors.push_back(flavor);
  }
  // throwaway addresses
  double csvWgtHF, csvWgtLF, csvWgtCF, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF;
  if (analysis->useCMVA) {
    gt->sf_cmvaweight = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,0, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF);
    gt->sf_cmva_LFup = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,9, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_LFdown = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,10, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_HFup = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,11, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_HFdown = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,12, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_HFStats1up = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,13, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_HFStats1down = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,14, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_LFStats1up = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,17, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_LFStats1down = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,18, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_HFStats2up = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,15, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_HFStats2down = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,16, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_LFStats2up = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,19, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_LFStats2down = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,20, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_CErr1up = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,21, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_CErr1down = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,22, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_CErr2up = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,23, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
    gt->sf_cmva_CErr2down = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors,24, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF)/gt->sf_cmvaweight;
  } else {
    gt->sf_csvweight = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,0, csvWgtHF, csvWgtLF, csvWgtCF);
    gt->sf_csv_LFup = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,9, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_LFdown = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,10, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_HFup = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,11, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_HFdown = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,12, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_HFStats1up = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,13, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_HFStats1down = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,14, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_LFStats1up = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,17, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_LFStats1down = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,18, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_HFStats2up = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,15, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_HFStats2down = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,16, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_LFStats2up = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,19, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_LFStats2down = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,20, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_CErr1up = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,21, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_CErr1down = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,22, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_CErr2up = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,23, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
    gt->sf_csv_CErr2down = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors,24, csvWgtHF, csvWgtLF, csvWgtCF)/gt->sf_csvweight;
  }
}
