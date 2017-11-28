#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define EGMSCALE 1

using namespace panda;
using namespace std;


float PandaAnalyzer::GetMSDCorr(Float_t puppipt, Float_t puppieta) 
{

  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;

  genCorr = puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta) <= 1.3 ){
    recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  }
  else{
    recoCorr = puppisd_corrRECO_for->Eval( puppipt );
  }
  totalWeight = genCorr * recoCorr;

  return totalWeight;
}

void PandaAnalyzer::FatjetBasics() 
{
  fj1=0;
  gt->nFatjet=0;
  int fatjet_counter=-1;
  for (auto& fj : *fatjets) {
    ++fatjet_counter;
    float pt = fj.pt();
    float rawpt = fj.rawPt;
    float eta = fj.eta();
    float mass = fj.m();
    float ptcut = 200;
    if (analysis->monoh)
      ptcut = 200;

    if (pt<ptcut || fabs(eta)>2.4 || !fj.monojet)
      continue;

    float phi = fj.phi();
    if (IsMatched(&matchLeps,FATJETMATCHDR2,eta,phi) || IsMatched(&matchPhos,FATJETMATCHDR2,eta,phi)) {
      continue;
    }

    gt->nFatjet++;
    if (gt->nFatjet==1) {
      fj1 = &fj;
      if (fatjet_counter==0)
        gt->fj1IsClean = 1;
      else
        gt->fj1IsClean = 0;
      gt->fj1Pt = pt;
      gt->fj1Eta = eta;
      gt->fj1Phi = phi;
      gt->fj1M = mass;
      gt->fj1MSD = fj.mSD;
      gt->fj1RawPt = rawpt;

      // do a bit of jet energy scaling
      if (analysis->varyJES) {
        double scaleUnc = (fj.ptCorrUp - gt->fj1Pt) / gt->fj1Pt; 
        gt->fj1PtScaleUp   = gt->fj1Pt  * (1 + 2*scaleUnc);
        gt->fj1PtScaleDown  = gt->fj1Pt  * (1 - 2*scaleUnc);
        gt->fj1MSDScaleUp  = gt->fj1MSD * (1 + 2*scaleUnc);
        gt->fj1MSDScaleDown = gt->fj1MSD * (1 - 2*scaleUnc);

        // do some jet energy smearing
        if (isData) {
          gt->fj1PtSmeared = gt->fj1Pt;
          gt->fj1PtSmearedUp = gt->fj1Pt;
          gt->fj1PtSmearedDown = gt->fj1Pt;
          gt->fj1MSDSmeared = gt->fj1MSD;
          gt->fj1MSDSmearedUp = gt->fj1MSD;
          gt->fj1MSDSmearedDown = gt->fj1MSD;
        } else {
          double smear=1, smearUp=1, smearDown=1;
          ak8JERReader->getStochasticSmear(pt,eta,event.rho,smear,smearUp,smearDown);

          gt->fj1PtSmeared = smear*gt->fj1Pt;
          gt->fj1PtSmearedUp = smearUp*gt->fj1Pt;
          gt->fj1PtSmearedDown = smearDown*gt->fj1Pt;

          gt->fj1MSDSmeared = smear*gt->fj1MSD;
          gt->fj1MSDSmearedUp = smearUp*gt->fj1MSD;
          gt->fj1MSDSmearedDown = smearDown*gt->fj1MSD;
        }

        // now have to do this mess with the subjets...
        TLorentzVector sjSum, sjSumUp, sjSumDown, sjSumSmear;
        for (unsigned int iSJ=0; iSJ!=fj.subjets.size(); ++iSJ) {
          auto& subjet = fj.subjets.objAt(iSJ);
          // now correct...
          double factor=1;
          if (fabs(subjet.eta())<5.191) {
            scaleReaderAK4->setJetPt(subjet.pt());
            scaleReaderAK4->setJetEta(subjet.eta());
            scaleReaderAK4->setJetPhi(subjet.phi());
            scaleReaderAK4->setJetE(subjet.e());
            scaleReaderAK4->setRho(event.rho);
            scaleReaderAK4->setJetA(0);
            scaleReaderAK4->setJetEMF(-99.0);
            factor = scaleReaderAK4->getCorrection();
          }
          TLorentzVector vCorr = factor * subjet.p4();
          sjSum += vCorr;
          double corr_pt = vCorr.Pt();

          // now vary
          uncReaderAK4->setJetEta(subjet.eta()); uncReaderAK4->setJetPt(corr_pt);
          double scaleUnc = uncReaderAK4->getUncertainty(true);
          sjSumUp += (1 + 2*scaleUnc) * vCorr;
          sjSumDown += (1 - 2*scaleUnc) * vCorr;

          // now smear...
          double smear=1, smearUp=1, smearDown=1;
          ak4JERReader->getStochasticSmear(corr_pt,subjet.eta(),event.rho,smear,smearUp,smearDown);
          sjSumSmear += smear * vCorr;
        }
        gt->fj1PtScaleUp_sj = gt->fj1Pt * (sjSumUp.Pt()/sjSum.Pt());
        gt->fj1PtScaleDown_sj = gt->fj1Pt * (sjSumDown.Pt()/sjSum.Pt());
        gt->fj1PtSmeared_sj = gt->fj1Pt * (sjSumSmear.Pt()/sjSum.Pt());
        gt->fj1MSDScaleUp_sj = gt->fj1MSD * (sjSumUp.Pt()/sjSum.Pt());
        gt->fj1MSDScaleDown_sj = gt->fj1MSD * (sjSumDown.Pt()/sjSum.Pt());
        gt->fj1MSDSmeared_sj = gt->fj1MSD * (sjSumSmear.Pt()/sjSum.Pt());
      }

      if (analysis->monoh) {
        // mSD correction
        float corrweight=1.;
        corrweight = GetMSDCorr(pt,eta);
        gt->fj1MSD_corr = corrweight*gt->fj1MSD;
      }

      // now we do substructure
      gt->fj1Tau32 = clean(fj.tau3/fj.tau2);
      gt->fj1Tau32SD = clean(fj.tau3SD/fj.tau2SD);
      gt->fj1Tau21 = clean(fj.tau2/fj.tau1);
      gt->fj1Tau21SD = clean(fj.tau2SD/fj.tau1SD);

      for (auto ibeta : ibetas) {
        for (auto N : Ns) {
          for (auto order : orders) {
            GeneralTree::ECFParams p;
            p.order = order; p.N = N; p.ibeta = ibeta;
            if (gt->fj1IsClean || true)
              gt->fj1ECFNs[p] = fj.get_ecf(order,N,ibeta);
            else
              gt->fj1ECFNs[p] = fj.get_ecf(order,N,ibeta);
          }
        }
      } //loop over betas
      gt->fj1HTTMass = fj.htt_mass;
      gt->fj1HTTFRec = fj.htt_frec;

      std::vector<panda::MicroJet const*> subjets;
      for (unsigned iS(0); iS != fj.subjets.size(); ++iS)
        subjets.push_back(&fj.subjets.objAt(iS));

      auto csvsort([](panda::MicroJet const* j1, panda::MicroJet const* j2)->bool {
          return j1->csv > j2->csv;
          });

      std::sort(subjets.begin(),subjets.end(),csvsort);
      if (subjets.size()>0) {
        gt->fj1MaxCSV = subjets.at(0)->csv;
        gt->fj1MinCSV = subjets.back()->csv;
        if (subjets.size()>1) {
          gt->fj1SubMaxCSV = subjets.at(1)->csv;
        }
      }

      if (analysis->monoh) {
        gt->fj1DoubleCSV = fj.double_sub;
        for (unsigned int iSJ=0; iSJ!=fj.subjets.size(); ++iSJ) {
          auto& subjet = fj.subjets.objAt(iSJ);
          gt->fj1sjPt[iSJ]=subjet.pt();
          gt->fj1sjEta[iSJ]=subjet.eta();
          gt->fj1sjPhi[iSJ]=subjet.phi();
          gt->fj1sjM[iSJ]=subjet.m();
          gt->fj1sjCSV[iSJ]=subjet.csv;
          gt->fj1sjQGL[iSJ]=subjet.qgl;
        }
      }
    }
  }
  tr->TriggerSubEvent("fatjet basics");
}

void PandaAnalyzer::FatjetRecluster() 
{
  if (fj1) {
    VPseudoJet particles = ConvertPFCands(event.pfCandidates,analysis->puppi_jets,0);
    fastjet::ClusterSequenceArea seq(particles,*jetDef,*areaDef);
    VPseudoJet allJets(seq.inclusive_jets(0.));
    fastjet::PseudoJet *pj1=0;
    double minDR2 = 999;
    for (auto &jet : allJets) {
      double dr2 = DeltaR2(jet.eta(),jet.phi_std(),fj1->eta(),fj1->phi());
      if (dr2<minDR2) {
        minDR2 = dr2;
        pj1 = &jet;
      }
    }
    if (pj1) {
      VPseudoJet constituents = fastjet::sorted_by_pt(pj1->constituents());

      gt->fj1NConst = constituents.size();
      double eTot=0, eTrunc=0;
      for (unsigned iC=0; iC!=gt->fj1NConst; ++iC) {
        double e = constituents.at(iC).E();
        eTot += e;
        if (iC<100)
          eTrunc += e;
      }
      gt->fj1EFrac100 = eTrunc/eTot;


      fastjet::PseudoJet sdJet = (*softDrop)(*pj1);
      VPseudoJet sdConstituents = fastjet::sorted_by_pt(sdJet.constituents());
      gt->fj1NSDConst = sdConstituents.size();
      eTot=0; eTrunc=0;
      for (unsigned iC=0; iC!=gt->fj1NSDConst; ++iC) {
        double e = sdConstituents.at(iC).E();
        eTot += e;
        if (iC<100)
          eTrunc += e;
      }
      gt->fj1SDEFrac100 = eTrunc/eTot;

    }
    tr->TriggerSubEvent("fatjet reclustering");
  }

}

void PandaAnalyzer::FatjetMatching() 
{
  // identify interesting gen particles for fatjet matching
  unsigned int pdgidTarget=0;
  if (!isData && analysis->processType>=kTT) {
    switch(analysis->processType) {
      case kTop:
      case kTT:
      case kSignal:
        pdgidTarget=6;
        break;
      case kV:
        pdgidTarget=24;
        break;
      case kH:
        pdgidTarget=25;
        break;
      default:
        // analysis->processType>=kTT means we should never get here
        PError("PandaAnalyzer::Run","Reached an unknown process type");
    }

    std::vector<int> targets;

    int nGen = event.genParticles.size();
    for (int iG=0; iG!=nGen; ++iG) {
      auto& part(event.genParticles.at(iG));
      int pdgid = part.pdgid;
      unsigned int abspdgid = abs(pdgid);
      if (abspdgid == pdgidTarget)
        targets.push_back(iG);
    } //looking for targets

    for (int iG : targets) {
      auto& part(event.genParticles.at(iG));

      // check there is no further copy:
      bool isLastCopy=true;
      for (int jG : targets) {
        if (event.genParticles.at(jG).parent.get() == &part) {
          isLastCopy=false;
          break;
        }
      }
      if (!isLastCopy)
        continue;

      // (a) check it is a hadronic decay and if so, (b) calculate the size
      if (analysis->processType==kTop||analysis->processType==kTT) {

        // first look for a W whose parent is the top at iG, or a W further down the chain
        panda::GenParticle const* lastW(0);
        for (int jG=0; jG!=nGen; ++jG) {
          GenParticle const& partW(event.genParticles.at(jG));
          if (TMath::Abs(partW.pdgid)==24 && partW.pdgid*part.pdgid>0) {
            // it's a W and has the same sign as the top
            if (!lastW && partW.parent.get() == &part) {
              lastW = &partW;
            } else if (lastW && partW.parent.get() == lastW) {
              lastW = &partW;
            }
          }
        } // looking for W
        if (!lastW) {// ???
          continue;
        }
        auto& partW(*lastW);

        // now look for b or W->qq
        int iB=-1, iQ1=-1, iQ2=-1;
        double size=0, sizeW=0;
        for (int jG=0; jG!=nGen; ++jG) {
          auto& partQ(event.genParticles.at(jG));
          int pdgidQ = partQ.pdgid;
          unsigned int abspdgidQ = TMath::Abs(pdgidQ);
          if (abspdgidQ>5)
            continue;
          if (abspdgidQ==5 && iB<0 && partQ.parent.get() == &part) {
            // only keep first copy
            iB = jG;
            size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),size);
          } else if (abspdgidQ<5 && partQ.parent.get() == &partW) {
            if (iQ1<0) {
              iQ1 = jG;
              size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),
                  size);
              sizeW = TMath::Max(DeltaR2(partW.eta(),partW.phi(),partQ.eta(),partQ.phi()),
                  sizeW);
            } else if (iQ2<0) {
              iQ2 = jG;
              size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),
                  size);
              sizeW = TMath::Max(DeltaR2(partW.eta(),partW.phi(),partQ.eta(),partQ.phi()),
                  sizeW);
            }
          }
          if (iB>=0 && iQ1>=0 && iQ2>=0)
            break;
        } // looking for quarks


        bool isHadronic = (iB>=0 && iQ1>=0 && iQ2>=0); // all 3 quarks were found
        if (isHadronic)
          genObjects[&part] = size;

        bool isHadronicW = (iQ1>=0 && iQ2>=0);
        if (isHadronicW)
          genObjects[&partW] = sizeW;

      } else { // these are W,Z,H - 2 prong decays

        int iQ1=-1, iQ2=-1;
        double size=0;
        for (int jG=0; jG!=nGen; ++jG) {
          auto& partQ(event.genParticles.at(jG));
          int pdgidQ = partQ.pdgid;
          unsigned int abspdgidQ = TMath::Abs(pdgidQ);
          if (abspdgidQ>5)
            continue;
          if (partQ.parent.get() == &part) {
            if (iQ1<0) {
              iQ1=jG;
              size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),
                  size);
            } else if (iQ2<0) {
              iQ2=jG;
              size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),
                  size);
            }
          }
          if (iQ1>=0 && iQ2>=0)
            break;
        } // looking for quarks

        bool isHadronic = (iQ1>=0 && iQ2>=0); // both quarks were found

        // add to collection
        if (isHadronic)
          genObjects[&part] = size;
      }

    } // loop over targets
  } // process is interesting

  tr->TriggerEvent("gen matching");

  if (!isData && gt->nFatjet>0) {
    // first see if jet is matched
    auto* matched = MatchToGen(fj1->eta(),fj1->phi(),1.5,pdgidTarget);
    if (matched!=NULL) {
      gt->fj1IsMatched = 1;
      gt->fj1GenPt = matched->pt();
      gt->fj1GenSize = genObjects[matched];
    } else {
      gt->fj1IsMatched = 0;
    }
    if (pdgidTarget==6) { // matched to top; try for W
      auto* matchedW = MatchToGen(fj1->eta(),fj1->phi(),1.5,24);
      if (matchedW!=NULL) {
        gt->fj1IsWMatched = 1;
        gt->fj1GenWPt = matchedW->pt();
        gt->fj1GenWSize = genObjects[matchedW];
      } else {
        gt->fj1IsWMatched = 0;
      }
    }

    bool found_b_from_g=false;
    int bs_inside_cone=0;
    int has_gluon_splitting=0;
    panda::GenParticle const* first_b_mo(0);
    // now get the highest pT gen particle inside the jet cone
    for (auto& gen : event.genParticles) {
      float pt = gen.pt();
      int pdgid = gen.pdgid;
      if (pt>(gt->fj1HighestPtGenPt)
          && DeltaR2(gen.eta(),gen.phi(),fj1->eta(),fj1->phi())<FATJETMATCHDR2) {
        gt->fj1HighestPtGenPt = pt;
        gt->fj1HighestPtGen = pdgid;
      }

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

      if (DeltaR2(gen.eta(),gen.phi(),fj1->eta(),fj1->phi())<FATJETMATCHDR2) {
        gt->fj1NHF++;
        if (apdgid==5) {
          if (gen.parent.isValid() && gen.parent->pdgid==21 && gen.parent->pt()>20) {
            if (!found_b_from_g) {
              found_b_from_g=true;
              first_b_mo=gen.parent.get();
              bs_inside_cone+=1;
            } else if (gen.parent.get()==first_b_mo) {
              bs_inside_cone+=1;
              has_gluon_splitting=1;
            } else {
              bs_inside_cone+=1;
            }
          } else {
            bs_inside_cone+=1;
          }
        }
      }
    }

    gt->fj1Nbs=bs_inside_cone;
    gt->fj1gbb=has_gluon_splitting;

    // now get the subjet btag SFs
    vector<btagcand> sj_btagcands;
    vector<double> sj_sf_cent, sj_sf_bUp, sj_sf_bDown, sj_sf_mUp, sj_sf_mDown;
    unsigned int nSJ = fj1->subjets.size();
    for (unsigned int iSJ=0; iSJ!=nSJ; ++iSJ) {
      auto& subjet = fj1->subjets.objAt(iSJ);
      int flavor=0;
      for (auto& gen : event.genParticles) {
        int apdgid = abs(gen.pdgid);
        if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
          continue;
        double dr2 = DeltaR2(subjet.eta(),subjet.phi(),gen.eta(),gen.phi());
        if (dr2<0.09) {
          if (apdgid==4 || apdgid==5) {
            flavor=apdgid;
            break;
          } else {
            flavor=0;
          }
        }
      } // finding the subjet flavor

      float pt = subjet.pt();
      float btagUncFactor = 1;
      float eta = subjet.eta();
      double eff(1),sf(1),sfUp(1),sfDown(1);
      unsigned int binpt = btagpt.bin(pt);
      unsigned int bineta = btageta.bin(fabs(eta));
      if (flavor==5) {
        eff = beff[bineta][binpt];
      } else if (flavor==4) {
        eff = ceff[bineta][binpt];
      } else {
        eff = lfeff[bineta][binpt];
      }
      CalcBJetSFs(bSubJetL,flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
      sj_btagcands.push_back(btagcand(iSJ,flavor,eff,sf,sfUp,sfDown));
      sj_sf_cent.push_back(sf);
      if (flavor>0) {
        sj_sf_bUp.push_back(sfUp); sj_sf_bDown.push_back(sfDown);
        sj_sf_mUp.push_back(sf); sj_sf_mDown.push_back(sf);
      } else {
        sj_sf_bUp.push_back(sf); sj_sf_bDown.push_back(sf);
        sj_sf_mUp.push_back(sfUp); sj_sf_mDown.push_back(sfDown);
      }

    } // loop over subjets

    EvalBTagSF(sj_btagcands,sj_sf_cent,GeneralTree::bCent,GeneralTree::bSubJet);
    EvalBTagSF(sj_btagcands,sj_sf_bUp,GeneralTree::bBUp,GeneralTree::bSubJet);
    EvalBTagSF(sj_btagcands,sj_sf_bDown,GeneralTree::bBDown,GeneralTree::bSubJet);
    EvalBTagSF(sj_btagcands,sj_sf_mUp,GeneralTree::bMUp,GeneralTree::bSubJet);
    EvalBTagSF(sj_btagcands,sj_sf_mDown,GeneralTree::bMDown,GeneralTree::bSubJet);

  }

  tr->TriggerEvent("fatjet gen-matching");
}

