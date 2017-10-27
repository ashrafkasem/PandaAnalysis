#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define EGMSCALE 1

using namespace panda;
using namespace std;


void PandaAnalyzer::TopPTReweight()
{
      if (processType != kTT)
        return;

      gt->genWPlusPt = -1; gt->genWMinusPt = -1;
      for (auto& gen : event.genParticles) {
        if (abs(gen.pdgid)!=24)
          continue;
        if (analysis->firstGen) {
          if (gen.parent.isValid() && gen.parent->pdgid==gen.pdgid)
            continue; // must be first copy
        }
        if (gen.pdgid>0) {
         gt->genWPlusPt = gen.pt();
         gt->genWPlusEta = gen.eta();
        } else {
         gt->genWMinusPt = gen.pt();
         gt->genWMinusEta = gen.eta();
        }
        if (analysis->firstGen) {
          if (gt->genWPlusPt>0 && gt->genWMinusPt>0)
            break;
        }
      }
      TLorentzVector vT,vTbar;
      float pt_t=0, pt_tbar=0;
      for (auto& gen : event.genParticles) {
        if (abs(gen.pdgid)!=6)
          continue;
        if (analysis->firstGen) {
          if (gen.parent.isValid() && gen.parent->pdgid==gen.pdgid)
            continue; // must be first copy
        }
        if (gen.pdgid>0) {
         pt_t = gen.pt();
         gt->genTopPt = gen.pt();
         gt->genTopEta = gen.eta();
         vT.SetPtEtaPhiM(gen.pt(),gen.eta(),gen.phi(),gen.m());
        } else {
         pt_tbar = gen.pt();
         gt->genAntiTopPt = gen.pt();
         gt->genAntiTopEta = gen.eta();
         vTbar.SetPtEtaPhiM(gen.pt(),gen.eta(),gen.phi(),gen.m());
        }
        if (analysis->firstGen) {
          if (pt_t>0 && pt_tbar>0)
            break;
        }
      }
      if (pt_t>0 && pt_tbar>0) {
        TLorentzVector vTT = vT+vTbar;
        gt->genTTPt = vTT.Pt(); gt->genTTEta = vTT.Eta();
        gt->sf_tt8TeV       = TMath::Sqrt(TMath::Exp(0.156-0.00137*TMath::Min((float)400.,pt_t)) *
                         TMath::Exp(0.156-0.00137*TMath::Min((float)400.,pt_tbar)));
        gt->sf_tt           = TMath::Sqrt(TMath::Exp(0.0615-0.0005*TMath::Min((float)400.,pt_t)) *
                         TMath::Exp(0.0615-0.0005*TMath::Min((float)400.,pt_tbar)));
        gt->sf_tt8TeV_ext   = TMath::Sqrt(TMath::Exp(0.156-0.00137*pt_t) *
                         TMath::Exp(0.156-0.00137*pt_tbar));
        gt->sf_tt_ext       = TMath::Sqrt(TMath::Exp(0.0615-0.0005*pt_t) *
                         TMath::Exp(0.0615-0.0005*pt_tbar));
        gt->sf_tt8TeV_bound = TMath::Sqrt(((pt_t>400) ? 1 : TMath::Exp(0.156-0.00137*pt_t)) *
                         ((pt_tbar>400) ? 1 : TMath::Exp(0.156-0.00137*pt_tbar)));
        gt->sf_tt_bound     = TMath::Sqrt(((pt_t>400) ? 1 : TMath::Exp(0.0615-0.0005*pt_t)) *
                         ((pt_tbar>400) ? 1 : TMath::Exp(0.0615-0.0005*pt_tbar)));
      }

      if (pt_t>0)
        gt->sf_qcdTT *= TTNLOToNNLO(pt_t);
      if (pt_tbar>0) 
        gt->sf_qcdTT *= TTNLOToNNLO(pt_tbar);
      gt->sf_qcdTT = TMath::Sqrt(gt->sf_qcdTT);

    tr->TriggerEvent("tt SFs");
}

void PandaAnalyzer::VJetsReweight() 
{
      // calculate the mjj 
      TLorentzVector vGenJet;
      if (analysis->vbf && event.ak4GenJets.size() > 0) {
        auto &gj = event.ak4GenJets.at(0);
        TLorentzVector v;
        v.SetPtEtaPhiM(gj.pt(), gj.eta(), gj.phi(), gj.m());
        gt->genJet1Pt = gj.pt(); gt->genJet1Eta = gj.eta();
        vGenJet += v;
        if (event.ak4GenJets.size() > 1) {
          gj = event.ak4GenJets.at(1);
          v.SetPtEtaPhiM(gj.pt(), gj.eta(), gj.phi(), gj.m());
          gt->genJet2Pt = gj.pt(); gt->genJet2Eta = gj.eta();
          vGenJet += v;
        }
      }
      gt->genMjj = vGenJet.M();

      bool found = processType!=kA 
                   && processType!=kZ 
                   && processType!=kW
                   && processType!=kZEWK 
                   && processType!=kWEWK;
      if (found)
        return;

      int target=24;
      if (processType==kZ || processType==kZEWK) target=23;
      if (processType==kA) target=22;

      for (auto& gen : event.genParticles) {
        if (found) break;
        int apdgid = abs(gen.pdgid);
        if (apdgid==target)     {
          bool foundChild = false;
          for (auto& child : event.genParticles) {
            if (abs(child.pdgid) != target)
              continue;
            if (child.parent.isValid() && child.parent.get() == &(gen)) {
              foundChild = true; 
              break;
            }
          }    
          if (foundChild)
            continue;
          if (processType==kZ) {
            gt->trueGenBosonPt = gen.pt();
            gt->genBosonMass = gen.m();
            gt->genBosonEta = gen.eta();
            gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
            gt->sf_qcdV = GetCorr(cZNLO,gt->genBosonPt);
            gt->sf_ewkV = GetCorr(cZEWK,gt->genBosonPt);
            if (analysis->vbf) {
              gt->sf_qcdV_VBF = GetCorr(cVBF_ZNLO,gt->genBosonPt,gt->genMjj);
              gt->sf_qcdV_VBF2l = GetCorr(cVBF_ZllNLO,gt->genBosonPt,gt->genMjj);
              gt->sf_qcdV_VBFTight = GetCorr(cVBFTight_ZNLO,gt->genBosonPt);
              gt->sf_qcdV_VBF2lTight = GetCorr(cVBFTight_ZllNLO,gt->genBosonPt);
            }
            found=true;
          } else if (processType==kW) {
            gt->trueGenBosonPt = gen.pt();
            gt->genBosonMass = gen.m();
            gt->genBosonEta = gen.eta();
            gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
            gt->sf_qcdV = GetCorr(cWNLO,gt->genBosonPt);
            gt->sf_ewkV = GetCorr(cWEWK,gt->genBosonPt);
            if (analysis->vbf) {
              gt->sf_qcdV_VBF = GetCorr(cVBF_WNLO,gt->genBosonPt,gt->genMjj);
              gt->sf_qcdV_VBFTight = GetCorr(cVBFTight_WNLO,gt->genBosonPt);
            }
            found=true;
          } else if (processType==kZEWK) {
            gt->trueGenBosonPt = gen.pt();
            gt->genBosonMass = gen.m();
            gt->genBosonEta = gen.eta();
            gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
            if (analysis->vbf) {
              gt->sf_qcdV_VBF = GetCorr(cVBF_EWKZ,gt->genBosonPt,gt->genMjj);
              gt->sf_qcdV_VBFTight = gt->sf_qcdV_VBF; // for consistency
            }
          } else if (processType==kWEWK) {
            gt->trueGenBosonPt = gen.pt();
            gt->genBosonMass = gen.m();
            gt->genBosonEta = gen.eta();
            gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
            if (analysis->vbf) {
              gt->sf_qcdV_VBF = GetCorr(cVBF_EWKW,gt->genBosonPt,gt->genMjj);
              gt->sf_qcdV_VBFTight = gt->sf_qcdV_VBF; // for consistency
            }
          } else if (processType==kA) {
            // take the highest pT
            if (gen.pt() > gt->trueGenBosonPt) {
              gt->trueGenBosonPt = gen.pt();
              gt->genBosonMass = gen.m();
              gt->genBosonEta = gen.eta();
              gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
              gt->sf_qcdV = GetCorr(cANLO,gt->genBosonPt);
              gt->sf_ewkV = GetCorr(cAEWK,gt->genBosonPt);
              gt->sf_qcdV2j = GetCorr(cANLO2j,gt->genBosonPt);
            }
          }
        } // target matches
      } // gen particle loop ends

      //now for the cases where we did not find a gen boson
      if (gt->genBosonPt < 0){

        TLorentzVector vpt(0,0,0,0);

        for (auto& part : event.genParticles) {
          int pdgid = part.pdgid;
          unsigned int abspdgid = abs(pdgid);

          if ((abspdgid == 11 || abspdgid == 13) &&
              (part.statusFlags == GenParticle::kIsPrompt || 
               part.statusFlags == GenParticle::kIsTauDecayProduct || 
               part.statusFlags == GenParticle::kIsPromptTauDecayProduct || 
               part.statusFlags == GenParticle::kIsDirectTauDecayProduct || 
               part.statusFlags == GenParticle::kIsDirectPromptTauDecayProduct )){

            //ideally you want to have dressed leptons (lepton + photon), 
            //but we have in any ways have a photon veto in the analysis
            if(IsMatched(&matchLeps,0.01,part.eta(),part.phi()))
              vpt += part.p4();
          }
          
          if ((abspdgid == 12 || abspdgid == 14 || abspdgid == 16) && part.finalState==1){
            vpt += part.p4();
          }
        }
        
        gt->genBosonPt = bound(vpt.Pt(),genBosonPtMin,genBosonPtMax);
        gt->trueGenBosonPt = vpt.Pt();
        gt->genBosonMass = vpt.M();
        gt->genBosonEta = vpt.Eta();

        if (processType==kZ) {
          gt->sf_qcdV = GetCorr(cZNLO,gt->genBosonPt);
          gt->sf_ewkV = GetCorr(cZEWK,gt->genBosonPt);
          if (analysis->vbf) {
            gt->sf_qcdV_VBF = GetCorr(cVBF_ZNLO,gt->genBosonPt,gt->genMjj);
            gt->sf_qcdV_VBF2l = GetCorr(cVBF_ZllNLO,gt->genBosonPt,gt->genMjj);
            gt->sf_qcdV_VBFTight = GetCorr(cVBFTight_ZNLO,gt->genBosonPt);
            gt->sf_qcdV_VBF2lTight = GetCorr(cVBFTight_ZllNLO,gt->genBosonPt);
          }
        } 
        else if (processType==kW) {
          gt->sf_qcdV = GetCorr(cWNLO,gt->genBosonPt);
          gt->sf_ewkV = GetCorr(cWEWK,gt->genBosonPt);
          if (analysis->vbf) {
            gt->sf_qcdV_VBF = GetCorr(cVBF_WNLO,gt->genBosonPt,gt->genMjj);
            gt->sf_qcdV_VBFTight = GetCorr(cVBFTight_WNLO,gt->genBosonPt);
          }
        } 
        else if (processType==kZEWK) {
          if (analysis->vbf) {
            gt->sf_qcdV_VBF = GetCorr(cVBF_EWKZ,gt->genBosonPt,gt->genMjj);
            gt->sf_qcdV_VBFTight = gt->sf_qcdV_VBF; // for consistency
          }
        } 
        else if (processType==kWEWK) {
          if (analysis->vbf) {
            gt->sf_qcdV_VBF = GetCorr(cVBF_EWKW,gt->genBosonPt,gt->genMjj);
            gt->sf_qcdV_VBFTight = gt->sf_qcdV_VBF; // for consistency
          }
        }  
      }

    tr->TriggerEvent("qcd/ewk SFs");
}


void PandaAnalyzer::SignalInfo()
{
      if (processType != kSignal)
        return;

      bool found=false, foundbar=false;
      TLorentzVector vMediator(0,0,0,0);
      for (auto& gen : event.genParticles) {
        if (found && foundbar)
          break;
        if (abs(gen.pdgid) != 18)
          continue;
        if (gen.parent.isValid() && gen.parent->pdgid == gen.pdgid)
          continue;
        if (gen.pdgid == 18 && !found) {
          found = true;
          vMediator += gen.p4();
        } else if (gen.pdgid == -18 && !foundbar) {
          foundbar = true;
          vMediator += gen.p4();
        }
      }
      if (found && foundbar) {
        gt->trueGenBosonPt = vMediator.Pt();
        gt->genBosonPt = bound(gt->trueGenBosonPt,175,1200);
      }
}

void PandaAnalyzer::QCDUncs()
{
      gt->pdfUp = 1 + event.genReweight.pdfDW;
      gt->pdfDown = 1 - event.genReweight.pdfDW;
      auto &genReweight = event.genReweight;
      for (unsigned iS=0; iS!=6; ++iS) {
        float s=1;
        switch (iS) {
          case 0:
            s = genReweight.r1f2DW; break;
          case 1:
            s = genReweight.r1f5DW; break;
          case 2:
            s = genReweight.r2f1DW; break;
          case 3:
            s = genReweight.r2f2DW; break;
          case 4:
            s = genReweight.r5f1DW; break;
          case 5:
            s = genReweight.r5f5DW; break;
          default:
            break;
        }
        gt->scale[iS] = s; 
        gt->scaleUp = max(float(gt->scaleUp),float(s));
        gt->scaleDown = min(float(gt->scaleDown),float(s));
      }
      tr->TriggerEvent("qcd uncertainties");
}

void PandaAnalyzer::SignalReweights()
{
      unsigned nW = wIDs.size();
      if (nW) {
        for (unsigned iW=0; iW!=nW; ++iW) {
          gt->signal_weights[wIDs[iW]] = event.genReweight.genParam[iW];
        }
      }
}
