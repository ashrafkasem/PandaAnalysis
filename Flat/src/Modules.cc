#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define EGMSCALE 1
#define FATJETMATCHDR2 2.25

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

float PandaAnalyzer::GetMSDCorr(Float_t puppipt, Float_t puppieta) 
{

  float genCorr   = 1.;
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

      if (gt->nLooseLep>1 && gt->looseLep1PdgId+gt->looseLep2PdgId==0) {
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
          gt->fj1PtScaleUp    = gt->fj1Pt  * (1 + 2*scaleUnc);
          gt->fj1PtScaleDown  = gt->fj1Pt  * (1 - 2*scaleUnc);
          gt->fj1MSDScaleUp   = gt->fj1MSD * (1 + 2*scaleUnc);
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

void PandaAnalyzer::JetBasics() 
{
    gt->barrelJet12Pt = 0;
    gt->barrelHT = 0;
    unsigned nBarrelJets = 0;
    panda::Jet *jet1=0, *jet2=0;
    jot1=0; jot2=0;
    jotUp1=0; jotUp2=0;
    jotDown1=0; jotDown2=0;
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


     if (jet.pt()>30) { // nominal jets
       cleanedJets.push_back(&jet);
       if (cleanedJets.size()<3) {
         bool isBad = GetCorr(cBadECALJets,jet.eta(),jet.phi()) > 0;
         if (isBad)
           gt->badECALFilter = 0;
       }
 
       float csv = (fabs(jet.eta())<2.5) ? jet.csv : -1;
       if (fabs(jet.eta())<2.4) {
         centralJets.push_back(&jet);
         if (centralJets.size()==1) {
           jet1 = &jet;
           gt->jet1Pt = jet.pt();
           gt->jet1Eta = jet.eta();
           gt->jet1Phi = jet.phi();
           gt->jet1CSV = csv;
           gt->jet1IsTight = jet.monojet ? 1 : 0;
         } else if (centralJets.size()==2) {
           jet2 = &jet;
           gt->jet2Pt = jet.pt();
           gt->jet2Eta = jet.eta();
           gt->jet2Phi = jet.phi();
           gt->jet2CSV = csv;
         }
       }

       vJet.SetPtEtaPhiM(jet.pt(),jet.eta(),jet.phi(),jet.m());
 
       if (analysis->vbf)
         JetVBFBasics(jet);

       if (analysis->monoh) 
         JetHbbBasics(jet);
 
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
         if (analysis->monoh) {
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

    if (analysis->vbf) {
      gt->nJot = cleanedJets.size();
      JetVBFSystem();
    }

    tr->TriggerEvent("jets");

}

void PandaAnalyzer::JetHbbBasics(panda::Jet& jet)
{
     float csv = (fabs(jet.eta())<2.5) ? jet.csv : -1;
     gt->jetPt[cleanedJets.size()-1]=jet.pt();
     gt->jetEta[cleanedJets.size()-1]=jet.eta();
     gt->jetPhi[cleanedJets.size()-1]=jet.phi();
     gt->jetE[cleanedJets.size()-1]=jet.m();
     gt->jetCSV[cleanedJets.size()-1]=csv;
     gt->jetQGL[cleanedJets.size()-1]=jet.qgl;
     tr->TriggerSubEvent("H->bb jet");
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
   // do jes variation OUTSIDE of pt>30 check
   if (jet.ptCorrUp>30) {
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
   }
   if (jet.ptCorrDown>30) {
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
      for (unsigned int i = 0;i<btaggedJets.size();i++){
        panda::Jet *jet_1 = btaggedJets.at(i);
        TLorentzVector hbbdaughter1;
        hbbdaughter1.SetPtEtaPhiM(jet_1->pt(),jet_1->eta(),jet_1->phi(),jet_1->m());
        for (unsigned int j = i+1;j<btaggedJets.size();j++){
          panda::Jet *jet_2 = btaggedJets.at(j);
          TLorentzVector hbbdaughter2;
          hbbdaughter2.SetPtEtaPhiM(jet_2->pt(),jet_2->eta(),jet_2->phi(),jet_2->m());
          TLorentzVector hbbsystem = hbbdaughter1 + hbbdaughter2;
          if (hbbsystem.Pt()>tmp_hbbpt){
            tmp_hbbpt = hbbsystem.Pt();
            tmp_hbbeta = hbbsystem.Eta();
            tmp_hbbphi = hbbsystem.Phi();
            tmp_hbbm = hbbsystem.M();
            tmp_hbbjtidx1 = btagindices.at(i);
            tmp_hbbjtidx2 = btagindices.at(j);
          }
        }
      }
      gt->hbbpt = tmp_hbbpt;
      gt->hbbeta = tmp_hbbeta;
      gt->hbbphi = tmp_hbbphi;
      gt->hbbm = tmp_hbbm;
      gt->hbbjtidx[0] = tmp_hbbjtidx1;
      gt->hbbjtidx[1] = tmp_hbbjtidx2;

      tr->TriggerEvent("monohiggs");
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

void PandaAnalyzer::FatjetMatching() 
{
    // identify interesting gen particles for fatjet matching
    unsigned int pdgidTarget=0;
    if (!isData && processType>=kTT) {
      switch(processType) {
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
          // processType>=kTT means we should never get here
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
        if (processType==kTop||processType==kTT) {

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
