#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define DEBUG 0
using namespace panda;
using namespace std;

PandaAnalyzer::PandaAnalyzer() {
}

PandaAnalyzer::~PandaAnalyzer() {
}

void PandaAnalyzer::ResetBranches() {
  genObjects.clear();
  matchPhos.clear();
  matchEles.clear();
  matchLeps.clear();
  gt->Reset();
}

void PandaAnalyzer::SetOutputFile(TString fOutName) {
  fOut = new TFile(fOutName,"RECREATE");
  tOut = new TTree("events","events");

  gt->Reset(); // to be extra safe and fill the map before setting addresses
  gt->WriteTree(tOut);

}

void PandaAnalyzer::Init(TTree *t, TTree *infotree)
{
  if (!t) return;
  tIn = t;
  t->SetBranchAddress("event",&event);
  t->SetBranchAddress("puppiCA15",&fatjets);
  t->SetBranchAddress("puppiAK4",&jets);
  t->SetBranchAddress("electron",&electrons);
  t->SetBranchAddress("muon",&muons);
  t->SetBranchAddress("tau",&taus);
  t->SetBranchAddress("photon",&photons);
  t->SetBranchAddress("pfmet",&pfmet);
  t->SetBranchAddress("puppimet",&puppimet);
  if (!isData) {
    t->SetBranchAddress("gen",&genparts);  
    TH1F *hDTotalMCWeight = new TH1F("hDTotalMCWeight","hDTotalMCWeight",4,-2,2);
    TString val("fabs(info.mcWeight)/(info.mcWeight)")
    infotree->Draw(val+">>hDTotalMCWeight",val);
    fOut->WriteTObject(hDTotalMCWeight,"hDTotalMCWeight");
  }

}

PGenParticle *PandaAnalyzer::MatchToGen(double eta, double phi, double radius) {
  PGenParticle *found=NULL;
  double r2 = radius*radius;

  unsigned int counter=0;
  for (map<PGenParticle*,float>::iterator iG=genObjects.begin();
        iG!=genObjects.end(); ++iG) {
    if (found!=NULL)
      break;
    if (DeltaR2(eta,phi,iG->first->eta,iG->first->phi)<r2) 
      found = iG->first;
  }

  return found;
}

void PandaAnalyzer::Terminate() {
  fOut->WriteTObject(tIn);
  fOut->WriteTObject(hDTotalMCWeight);
  fOut->Close();

  fEleTrigB->Close();
  fEleTrigE->Close();
  fPhoTrig->Close();
  fEleSF->Close();
  fMuSF->Close();
  fEleSFTight->Close();
  fMuSFTight->Close();
  fEleSFTrack->Close();
  fMuSFTrack->Close();
  fPU->Close();
  fKFactor->Close();
  fHFNew->Close();
  fLFNew->Close();

  delete btagCalib;
  delete hfReader;
  delete lfReader;
  delete hfUpReader;
  delete lfUpReader;
  delete hfDownReader;
  delete lfDownReader;

//  delete ak8jec;
//  delete ak8unc;
}

void PandaAnalyzer::SetDataDir(const char *s) {
  TString dirPath(s);

  fEleTrigB    = new TFile(dirPath+"/trigger_eff/ele_trig_lowpt_rebinned.root");
  fEleTrigE    = new TFile(dirPath+"/trigger_eff/ele_trig_lowpt_rebinned.root");
  fPhoTrig     = new TFile(dirPath+"/trigger_eff/pho_trig.root");
  fMetTrig     = new TFile(dirPath+"/trigger_eff/met_trig.root");
  fEleTrigLow  = new TFile(dirPath+"/trigger_eff/ele_trig_lowpt.root");

  fEleSF        = new TFile(dirPath+"/scaleFactor_electron_vetoid_12p9.root");
  fEleSFTight   = new TFile(dirPath+"/scaleFactor_electron_tightid_12p9.root");
  fEleSFTrack   = new TFile(dirPath+"/scaleFactor_electron_track.root");

  fMuSF         = new TFile(dirPath+"/scaleFactor_muon_looseid_12p9.root");
  fMuSFTight    = new TFile(dirPath+"/scaleFactor_muon_tightid_12p9.root");
  fMuSFTrack    = new TFile(dirPath+"/scaleFactor_muon_track.root");

  fPU      = new TFile(dirPath+"/puWeight_13invfb.root");

  fKFactor = new TFile(dirPath+"/kfactors.root");

  hEleTrigB = (TH1D*) fEleTrigB->Get("h_num");
  hEleTrigE = (TH1D*) fEleTrigE->Get("h_num_endcap");

  hPhoTrig = (TH1D*) fPhoTrig->Get("h_num");
  hMetTrig = (TH1D*) fMetTrig->Get("numer");
  hEleTrigLow = (TH2D*) fEleTrigLow->Get("hEffEtaPt");

  hEleVeto  = (TH2D*) fEleSF->Get("scaleFactor_electron_vetoid_RooCMSShape");
  hEleTight = (TH2D*) fEleSFTight->Get("scaleFactor_electron_tightid_RooCMSShape");

  hMuLoose = (TH2D*) fMuSF->Get("scaleFactor_muon_looseid_RooCMSShape");
  hMuTight = (TH2D*) fMuSFTight->Get("scaleFactor_muon_tightid_RooCMSShape");

  hMuTrack = (TH1D*) fMuSFTrack->Get("htrack2");
  hEleTrack = (TH2D*) fEleSFTrack->Get("EGamma_SF2D");

  hPUWeight = (TH1D*)fPU->Get("hPU");

  hZNLO = (TH1D*)fKFactor->Get("ZJets_012j_NLO/nominal");
  hWNLO = (TH1D*)fKFactor->Get("WJets_012j_NLO/nominal");
  hANLO = (TH1D*)fKFactor->Get("GJets_1j_NLO/nominal_G");

  hZLO  = (TH1D*)fKFactor->Get("ZJets_LO/inv_pt");
  hWLO  = (TH1D*)fKFactor->Get("WJets_LO/inv_pt");
  hALO  = (TH1D*)fKFactor->Get("GJets_LO/inv_pt_G");
 
  hZEWK = (TH1D*)fKFactor->Get("EWKcorr/Z");
  hWEWK = (TH1D*)fKFactor->Get("EWKcorr/W");
  hAEWK = (TH1D*)fKFactor->Get("EWKcorr/photon");

  hZEWK->Divide(hZNLO);   hWEWK->Divide(hWNLO);   hAEWK->Divide(hANLO);
  hZNLO->Divide(hZLO);    hWNLO->Divide(hWLO);    hANLO->Divide(hALO);

  btagCalib = new BTagCalibration("csvv2",(dirPath+"/CSVv2_ichep.csv").Data());
  hfReader = new BTagCalibrationReader(btagCalib,BTagEntry::OP_LOOSE,"comb","central");
  lfReader = new BTagCalibrationReader(btagCalib,BTagEntry::OP_LOOSE,"incl","central");
  hfUpReader = new BTagCalibrationReader(btagCalib,BTagEntry::OP_LOOSE,"comb","up");
  lfUpReader = new BTagCalibrationReader(btagCalib,BTagEntry::OP_LOOSE,"incl","up");
  hfDownReader = new BTagCalibrationReader(btagCalib,BTagEntry::OP_LOOSE,"comb","down");
  lfDownReader = new BTagCalibrationReader(btagCalib,BTagEntry::OP_LOOSE,"incl","down");

  sj_btagCalib = new BTagCalibration("csvv2",(dirPath+"/subjet_CSVv2_ichep.csv").Data());
  sj_hfReader = new BTagCalibrationReader(sj_btagCalib,BTagEntry::OP_LOOSE,"lt","central");
  sj_lfReader = new BTagCalibrationReader(sj_btagCalib,BTagEntry::OP_LOOSE,"incl","central");
  sj_hfUpReader = new BTagCalibrationReader(sj_btagCalib,BTagEntry::OP_LOOSE,"lt","up");
  sj_lfUpReader = new BTagCalibrationReader(sj_btagCalib,BTagEntry::OP_LOOSE,"incl","up");
  sj_hfDownReader = new BTagCalibrationReader(sj_btagCalib,BTagEntry::OP_LOOSE,"lt","down");
  sj_lfDownReader = new BTagCalibrationReader(sj_btagCalib,BTagEntry::OP_LOOSE,"incl","down");

//  ak8jec = new JetCorrectorParameters((dirPath+"/Spring16_25nsV6_MC_Uncertainty_AK8PFPuppi.txt").Data());
//  ak8unc = new JetCorrectionUncertainty(*ak8jec);
}

bool PandaAnalyzer::PassEventFilters() {
  return true;
}

bool PandaAnalyzer::PassGoodLumis() {
  return true;
}

bool PandaAnalyzer::PassPreselection() {
  if (preselBits==0)
    return true;
  bool isGood=false;
  if (preselBits & kMonotopCA15) {
    if (gt->nFatjet==1 && gt->fj1Pt>250) { 
      if ( (gt->puppimet>200 || gt->UZmag>200 || gt->UWmag>200 || gt->UAmag>200) ||
            (gt->pfmet>200 || gt->pfUZmag>200 || gt->pfUWmag>200 || gt->pfUAmag>200) ) {
            isGood = true;
      }
    }
  }
  // if (preselBits & kVBF) {
  //   if (nSelectedJet>1 && jet1Pt>40 && jet1IsTight==1) {
  //     if ( (met>180 || pfUZmag>180 || pfUWmag>180 || pfUAmag>180) ||
  //           (puppimet>180 || UZmag>180 || UWmag>180 || UAmag>180) ) {
  //       isGood = true;
  //     }
  //   }
  // }
  
  return isGood;
}

// run
void PandaAnalyzer::Run() {

  // INITIALIZE --------------------------------------------------------------------------
  
  unsigned int nEvents = tIn->GetEntries();
  unsigned int nZero = 0;
  if (lastEvent>=0 && lastEvent<(int)nEvents)
    nEvents = lastEvent;
  if (firstEvent>=0)
    nZero = firstEvent;

  if (!fOut || !tIn) {
    PError("PandaAnalyzer::Run","NOT SETUP CORRECTLY");
    exit(1);
  }

  // get bounds
  float sf_puMin = 0, sf_puMax=999;
  float sf_eleEtaMax=9999, sf_muEtaMax=9999;
  float sf_elePtMin=0, sf_elePtMax=9999;
  float sf_muPtMin=0, sf_muPtMax=9999;
  float genBosonPtMin=150, genBosonPtMax=1000;
  if (!isData) {
    sf_puMin = hPUWeight->GetBinCenter(1);
    sf_puMax = hPUWeight->GetBinCenter(hPUWeight->GetNbinsX());

    const TAxis *sf_eleEta = hEleTight->GetXaxis();
    sf_eleEtaMax = sf_eleEta->GetBinCenter(sf_eleEta->GetNbins());

    const TAxis *sf_elePt = hEleTight->GetYaxis();
    sf_elePtMin = sf_elePt->GetBinCenter(1); sf_elePtMax = sf_elePt->GetBinCenter(sf_elePt->GetNbins());

    const TAxis *sf_muEta = hMuTight->GetXaxis();
    sf_muEtaMax = sf_muEta->GetBinCenter(sf_muEta->GetNbins());

    const TAxis *sf_muPt = hMuTight->GetYaxis();
    sf_muPtMin = sf_muPt->GetBinCenter(1); sf_muPtMax = sf_muPt->GetBinCenter(sf_muPt->GetNbins());

    genBosonPtMin = hZNLO->GetBinCenter(1); genBosonPtMax = hZNLO->GetBinCenter(hZNLO->GetNbinsX());
  }

  // these are bins of b-tagging eff in pT
  std::vector<double> vbtagpt {50,70,100,140,200,300,670};
  std::vector<double> beff  {0.60592 , 0.634613, 0.645663, 0.646712, 0.649283, 0.621973, 0.554093};
  std::vector<double> ceff  {0.122067, 0.119659, 0.123723, 0.132141, 0.143654, 0.143127, 0.133581};
  std::vector<double> lfeff {0.014992, 0.012208, 0.011654, 0.011675, 0.015165, 0.016569, 0.020099};
  Binner btagpt(vbtagpt);

  std::vector<double> vnewbtagpt {30,40,60,100,160};
  std::vector<double> vnewbtageta {0.8,1.6,2.41};
  Binner newbtagpt(vnewbtagpt);
  Binner newbtageta(vnewbtageta);

  // these are triggers. at some point these ought to be read from the file
  std::vector<unsigned int> metTriggers {0,1,2,3,4,5,6};
  std::vector<unsigned int> eleTriggers {14,15,16,17,18};
  std::vector<unsigned int> phoTriggers {20,21,24};
  std::vector<unsigned int> muTriggers  {7,8,9,10,12,12,13};

  // set up reporters
  unsigned int iE=0;
  ProgressReporter pr("PandaAnalyzer::Run",&iE,&nEvents,10);
  TimeReporter tr("PandaAnalyzer::Run",DEBUG);

  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    tr.Start();
    pr.Report();
    ResetBranches();
    tIn->GetEntry(iE);
    tr.TriggerEvent("GetEntry");

    // event info
    gt->mcWeight = (event->mcWeight>0) ? 1 : -1;
    gt->runNumber = event->runNumber;
    gt->lumiNumber = event->lumiNumber;
    gt->eventNumber = event->eventNumber;
    gt->npv = event->npv;
    gt->metFilter = (event->metfilters->at(0)) ? 1 : 0;
    if (!isData) 
      gt->sf_pu = getVal(hPUWeight,bound(npv,sf_puMin,sf_puMax));

    tr.TriggerEvent("initialize");

    // identify interesting gen particles for fatjet matching
    if (processType>kTT) {
      std::vector<int> targets;

      int nGen = genparts->size();
      for (int iG=0; iG!=nGen; ++iG) {
        PGenParticle *part = genparts->at(iG);
        int pdgid = part->pdgid;
        unsigned int abspdgid = TMath::Abs(pdgid);
        bool good= (  (processType==kTop && abspdgid==6) ||
                      (processType==kV && (abspdgid==23 || abspdgid==24)) ||
                      (processType==kH && (abspdgid==25))
                    );
        if (good)
          targets.push_back(iG);
      } //looking for targets

      for (int iG : targets) {
        PGenParticle *part = genparts->at(iG);

        // check there is no further copy:
        bool isLastCopy=true;
        for (int jG : targets) {
          if (genparts->at(jG)->parent==iG) {
            isLastCopy=false;
            break;
          }
        }
        if (!isLastCopy)
          continue;
        
        // (a) check it is a hadronic decay and if so, (b) calculate the size
        if (processType==kTop) {
          
          // first look for a W whose parent is the top at iG, or a W further down the chain
          int iW=-1;
          for (int jG=0; jG!=nGen; ++jG) {
            PGenParticle *partW = genparts->at(jG);
            if (TMath::Abs(partW->pdgid)==24 && partW->pdgid*part->pdgid>0) {
              // it's a W and has the same sign as the top
              if (iW<0 && partW->parent==iG) {
                iW=jG;
              } else if (iW>=0 && partW->parent==iW) {
                iW=jG;
              }
            }
          } // looking for W

          // now look for b or W->qq
          int iB=-1, iQ1=-1, iQ2=-1;
          if (iW<0) // ???
            continue;
          double size=0;
          for (int jG=0; jG!=nGen; ++jG) {
            PGenParticle *partQ = genparts->at(jG);
            int pdgidQ = partQ->pdgid;
            unsigned int abspdgidQ = TMath::Abs(pdgidQ);
            if (abspdgidQ>5)
              continue;
            if (abspdgidQ==5 && iB<0 && partQ->parent==iG) {
              // only keep first copy
              iB = jG;
              size = TMath::Max(DeltaR2(part->eta,part->phi,partQ->eta,partQ->phi),size);
            } else if (abspdgidQ<5 && partQ->parent==iW) {
              if (iQ1<0) {
                iQ1 = jG;
                size = TMath::Max(DeltaR2(part->eta,part->phi,partQ->eta,partQ->phi),size);
              } else if (iQ2<0) {
                iQ2 = jG;
                size = TMath::Max(DeltaR2(part->eta,part->phi,partQ->eta,partQ->phi),size);
              }
            }
            if (iB>=0 && iQ1>=0 && iQ2>=0)
              break;
          } // looking for quarks

          bool isHadronic = (iB>=0 && iQ1>=0 && iQ2>=0); // all 3 quarks were found

          // add to collection
          if (isHadronic)
            genObjects[part] = size;

        } else { // these are W,Z,H - 2 prong decays

          int iQ1=-1, iQ2=-1;
          double size=0;
          for (int jG=0; jG!=nGen; ++jG) {
            PGenParticle *partQ = genparts->at(jG);
            int pdgidQ = partQ->pdgid;
            unsigned int abspdgidQ = TMath::Abs(pdgidQ);
            if (abspdgidQ>5)
              continue; 
            if (partQ->parent==iG) {
              if (iQ1<0) {
                iQ1=jG;
                size = TMath::Max(DeltaR2(part->eta,part->phi,partQ->eta,partQ->phi),size);
              } else if (iQ2<0) {
                iQ2=jG;
                size = TMath::Max(DeltaR2(part->eta,part->phi,partQ->eta,partQ->phi),size);
              }
            }
            if (iQ1>=0 && iQ2>=0)
              break;
          } // looking for quarks

          bool isHadronic = (iQ1>=0 && iQ2>=0); // both quarks were found

          // add to coll0ection
          if (isHadronic)
            genObjects[part] = size;
        }

      } // loop over targets
    } // process is interesting

    tr.TriggerEvent("gen matching");

    // met
    gt->pfmet = pfmet->pt;
    gt->pfmetphi = pfmet->phi;
    gt->puppimet = puppimet->pt;
    gt->puppimetphi = puppimet->phi;
    TLorentzVector vPFMET, vPuppiMET;
    vPFMET.SetPtEtaPhiM(gt->pfmet,0,gt->pfmetphi,0); 
    vPuppiMET.SetPtEtaPhiM(gt->puppimet,0,gt->puppimetphi,0); 

    tr.TriggerEvent("met");

    //electrons
    std::vector<PObject*> looseLeps, tightLeps;
    for (PElectron *ele : electrons) {
      float pt = ele->pt; float eta = ele->eta; float aeta = fabs(eta);
      if (pt<10 || aeta>2.5 || (aeta>1.4442 && aeta<1.566))
        continue;
      if ((ele->id&PElectron::kVeto)==0)
        continue;
      if (!ElectronIsolation(pt,eta,ele->iso,PElectron::kVeto))
        continue;
      looseLeps.push_back(ele);
      gt->nLooseElectron++;
    } 

    // muons
    for (PMuon *mu : muons) {
      float pt = mu->pt; float eta = mu->eta; float aeta = fabs(eta);
      if (pt<10 || aeta>2.4)
        continue;
      if ((mu->id&PMuon::kLoose)==0)
        continue;
      if (!MuonIsolation(pt,eta,mu->iso,PMuon::kLoose))
        continue;
      looseLeps.push_back(mu);
      gt->nLooseMuon++;
    } 

    // now consider all leptons
    gt->nLooseLep = looseLeps.size();
    std::partial_sort(looseLeps.begin(),looseLeps.begin()+3,looseLeps.end(),SortPObjects);
    int lep_counter=1;
    for (PObject *lep : looseLeps) {
      if (lep_counter==1) {
        gt->looseLep1Pt = lep->pt;
        gt->looseLep1Eta = lep->eta;
        gt->looseLep1Phi = lep->phi;
      } else if (lep_counter==2) {
        gt->looseLep2Pt = lep->pt;
        gt->looseLep2Eta = lep->eta;
        gt->looseLep2Phi = lep->phi;
      } else {
        break;
      }
      // now specialize lepton types
      PMuon *mu = dyanmic_cast<PMuon*>(lep);
      if (mu!=NULL) {
        bool isTight = ( (mu->id&PMuon::kTight)!=0 &&
                          MuonIsolation(mu->pt,mu->eta,mu->iso,PMuon::kTight) &&
                          mu->pt>20 && fabs(mu->eta)<2.4 );
        if (lep_counter==1) {
          gt->looseLep1PdgId = mu->q*13;
          if (isTight) {
            gt->nTightMuon++;
            gt->looseLep1IsTight = 1;
            matchLeps.push_back(lep);
          }
        } else if (lep_counter==2) {
          gt->looseLep2PdgId = mu->q*13;
          if (isTight) {
            gt->nTightMuon++;
            gt->looseLep2IsTight = 1;
          }
          if (isTight || gt->looseLep1IsTight)
            matchLeps.push_back(lep);
        }
      } else {
        PElectron *ele = dyanmic_cast<PElectron*>(lep);
        bool isTight = ( (ele->id&PElectron::kTight)!=0 &&
                          ElectronIsolation(ele->pt,ele->eta,ele->iso,PElectron::kTight) &&
                          ele->pt>40 && fabs(ele->eta)<2.5 );
        if (lep_counter==1) {
          gt->looseLep1PdgId = mu->q*11;
          if (isTight) {
            gt->nTightMuon++;
            gt->looseLep1IsTight = 1;
            matchLeps.push_back(lep);
            matchEles.push_back(lep);
          }
        } else if (lep_counter==2) {
          gt->looseLep2PdgId = mu->q*11;
          if (isTight) {
            gt->nTightMuon++;
            gt->looseLep2IsTight = 1;
          }
          if (isTight || gt->looseLep1IsTight) {
            matchLeps.push_back(lep);
            matchEles.push_back(lep);
          }
        }
      }
      ++lep_counter;
    }
    gt->nTightLep = gt->nTightElectron + gt->nTightMuon;
    if (gt->nLooseLep>1 && gt->looseLep1PdgId+gt->looseLep2PdgId=0) {
      TLorentzVector v1,v2;
      PObject *lep1=looseLeps[0], *lep2=looseLeps[1];
      v1.SetPtEtaPhiM(lep1->pt,lep1->eta,lep1->phi,lep1->m);
      v2.SetPtEtaPhiM(lep2->pt,lep2->eta,lep2->phi,lep2->m);
      gt->diLepMass = (v1+v2).M();
    } else {
      gt->diLepMass = -1;
    }

    tr.TriggerEvent("leptons") 

    std::vector<panda::PPhoton*> loosePhos;
    for (PPhoton *pho : photons) {
      if ((pho->id&PPhoton::kLoose)==0)
        continue;
      float pt = pho->pt;
      if (pt<1) continue;
      float eta = pho->eta, phi = pho->phi;
      if (pt<15 || fabs(eta)>2.5)
        continue;
      if (IsMatched(&matchEles,0,0.16,eta,phi))
        continue;
      loosePhos.push_back(pho);
      gt->nLoosePhoton++;
      if (gt->nLoosePhoton==1) {
        gt->loosePho1Pt = pt;
        gt->loosePho1Eta = eta;
        gt->loosePho1Phi = phi;
      }
      if ( (pho->id&PPhoton::kTight)!=0 &&
            pt>175 && fabs(eta)<1.4442 ) {
        if (gt->nLoosePhoton==1)
          gt->loosePho1IsTight=1;
        gt->nTightPhoton++;
        matchPhos.push_back(pho);
      }
    }

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

    tr.TriggerEvent("photons");

    //HERE I AM - do recol next

    // these are ak4 jets stored in the tree
    for (auto *anajet : anajets) {
      JetWriter *outjet = anajet->outjet;
      VJet *injets = anajet->injets;

      int nJets = injets->size();
      for (int iJ=0; iJ!=nJets; ++iJ) {
        outjet->reset();
        PJet *pjet = injets->at(iJ);
        outjet->read(pjet);
        outjet->idx = iJ;
        anajet->outtree->Fill();
      } // loop over jets
    } // loop over jet collections

    tr.TriggerEvent("jets");

    // these are fat jets stored in the tree
    // much more complicated!
    for (auto *anafatjet : anafatjets) {
      FatJetWriter *outjet = anafatjet->outjet;
      VFatJet *injets = anafatjet->injets;

      int nJets = injets->size();
      for (int iJ=0; iJ!=nJets; ++iJ) {
        outjet->reset();
        if (maxJets>=0 && iJ==maxJets)
          break;
        PFatJet *pfatjet = injets->at(iJ);
        if (pfatjet->pt<minFatJetPt)
          continue;
        outjet->read(pfatjet);
        outjet->idx = iJ;

        PGenParticle *matched = Match(pfatjet->eta,pfatjet->phi,anafatjet->radius);
        if (matched!=NULL) {
          outjet->matched = 1;
          outjet->genpt = matched->pt;
          outjet->gensize = genObjects[matched];
        } else { 
          outjet->matched = 0; 
        }

        tr.TriggerSubEvent("gen matching");

        /////// fastjet ////////
        VPseudoJet vpj = ConvertFatJet(pfatjet,anafatjet->pfcands,0.1);

        fastjet::ClusterSequenceArea seqCA(vpj, *(anafatjet->jetDefCA), *areaDef);
        fastjet::ClusterSequenceArea seqAK(vpj, *(anafatjet->jetDefAK), *areaDef);

        VPseudoJet alljetsCA(seqCA.inclusive_jets(0.));
        fastjet::PseudoJet *leadingJetCA=0; 
        for (auto &jet : alljetsCA) {
          if (!leadingJetCA || jet.perp2()>leadingJetCA->perp2())
            leadingJetCA = &jet;
        }

        VPseudoJet alljetsAK(seqAK.inclusive_jets(0.));
        fastjet::PseudoJet *leadingJetAK=0;
        for (auto &jet : alljetsAK) {
          if (!leadingJetAK || jet.perp2()>leadingJetAK->perp2())
            leadingJetAK = &jet;
        }
        
        tr.TriggerSubEvent("clustering");

        if (leadingJetCA!=NULL || leadingJetAK!=NULL) {
          fastjet::PseudoJet sdJetCA = (*anafatjet->sd)(*leadingJetCA);
          fastjet::PseudoJet sdJetAK = (*anafatjet->sd)(*leadingJetAK);
          
          VPseudoJet sdsubjets = fastjet::sorted_by_pt(sdJetCA.exclusive_subjets_up_to(3));

          // get the constituents and sort them
          VPseudoJet sdConstituentsCA = sdJetCA.constituents();
          std::sort(sdConstituentsCA.begin(),sdConstituentsCA.end(),orderPseudoJet);

          VPseudoJet sdConstituentsAK = sdJetAK.constituents();
          std::sort(sdConstituentsAK.begin(),sdConstituentsAK.end(),orderPseudoJet);

          tr.TriggerSubEvent("soft drop");

          /////////// let's calculate ECFs! ///////////

          // filter the constituents
          int nFilter;
          nFilter = TMath::Min(100,(int)sdConstituentsCA.size());
          VPseudoJet sdConstituentsCAFiltered(sdConstituentsCA.begin(),sdConstituentsCA.begin()+nFilter);
          
          if (doECF) {
            for (auto beta : betas) {
              calcECFN(beta,sdConstituentsCAFiltered,ecfnmanager);
              for (auto N : Ns) {
                for (auto o : orders) {
                  outjet->ecfns["ecfN_"+makeECFString(o,N,beta)] = ecfnmanager->ecfns[TString::Format("%i_%i",N,o)];
                }
              }
            }
            for (auto beta : betas) {
              calcECFN(beta,sdConstituentsCAFiltered,ecfnmanager,false);
              for (auto N : Ns) {
                for (auto o : orders) {
                  outjet->ecfns["maxecfN_"+makeECFString(o,N,beta)] = ecfnmanager->ecfns[TString::Format("%i_%i",N,o)];
                }
              }
            }

            tr.TriggerSubEvent("ecfns");

            // now we calculate ECFs for the subjets
            unsigned int nS = sdsubjets.size();
            for (auto beta : betas) {
              if (beta<2.0)
                continue; // speed things up for now
              for (unsigned int iS=0; iS!=nS; ++iS) {
                VPseudoJet subconstituents = sdsubjets[iS].constituents();
                nFilter = TMath::Min(80,(int)subconstituents.size());
                std::sort(subconstituents.begin(),subconstituents.end(),orderPseudoJet);
                VPseudoJet subconstituentsFiltered(subconstituents.begin(),subconstituents.begin()+nFilter);

                calcECFN(beta,subconstituentsFiltered,subecfnmanager);
                outjet->subecfns["min_secfN_"+makeECFString(1,3,beta)] = TMath::Min(
                    (double)subecfnmanager->ecfns[TString::Format("%i_%i",3,1)],
                    (double)outjet->subecfns["min_secfN_"+makeECFString(1,3,beta)]);
                outjet->subecfns["min_secfN_"+makeECFString(2,3,beta)] = TMath::Min(
                    (double)subecfnmanager->ecfns[TString::Format("%i_%i",3,2)],
                    (double)outjet->subecfns["min_secfN_"+makeECFString(2,3,beta)]);
                outjet->subecfns["min_secfN_"+makeECFString(3,3,beta)] = TMath::Min(
                    (double)subecfnmanager->ecfns[TString::Format("%i_%i",3,3)],
                    (double)outjet->subecfns["min_secfN_"+makeECFString(3,3,beta)]);
                outjet->subecfns["sum_secfN_"+makeECFString(1,3,beta)] += subecfnmanager->ecfns[TString::Format("%i_%i",3,1)];
                outjet->subecfns["sum_secfN_"+makeECFString(2,3,beta)] += subecfnmanager->ecfns[TString::Format("%i_%i",3,2)];
                outjet->subecfns["sum_secfN_"+makeECFString(3,3,beta)] += subecfnmanager->ecfns[TString::Format("%i_%i",3,3)];
              }
              outjet->subecfns["avg_secfN_"+makeECFString(1,3,beta)] = outjet->subecfns["sum_secfN_"+makeECFString(1,3,beta)]/nS;
              outjet->subecfns["avg_secfN_"+makeECFString(2,3,beta)] = outjet->subecfns["sum_secfN_"+makeECFString(2,3,beta)]/nS;
              outjet->subecfns["avg_secfN_"+makeECFString(3,3,beta)] = outjet->subecfns["sum_secfN_"+makeECFString(3,3,beta)]/nS;
            }

            tr.TriggerSubEvent("subecfns");

          }

          //////////// now let's do groomed tauN! /////////////
          double tau3 = anafatjet->tau->getTau(3,sdConstituentsCA);
          double tau2 = anafatjet->tau->getTau(2,sdConstituentsCA);
          double tau1 = anafatjet->tau->getTau(1,sdConstituentsCA);
          outjet->tau32SD = clean(tau3/tau2);
          outjet->tau21SD = clean(tau2/tau1);

          tr.TriggerSubEvent("tauSD");

          //////////// Q-jet quantities ////////////////
          if (doQjets) {
            std::vector<qjetwrapper> q_jets = getQjets(vpj,qplugin,qdef,qcounter++,15,anafatjet->tau);
            //std::vector<qjetwrapper> q_jets = getQjets(vpj,qplugin,qdef,qcounter++,10);

            JetQuantity getmass = [](qjetwrapper w) { return w.jet.m(); };
            outjet->qmass = clean(qVolQuantity(q_jets,getmass));

            JetQuantity getpt= [](qjetwrapper w) { return w.jet.pt(); };
            outjet->qpt = clean(qVolQuantity(q_jets,getpt));

            JetQuantity gettau32 = [](qjetwrapper w) { return w.tau32; };
            outjet->qtau32 = clean(qVolQuantity(q_jets,gettau32));

            JetQuantity gettau21 = [](qjetwrapper w) { return w.tau21; };
            outjet->qtau21 = clean(qVolQuantity(q_jets,gettau21));

            tr.TriggerSubEvent("qvol");
          }

          //////////// heat map! ///////////////
          if (doHeatMap) {
            outjet->hmap = HeatMap(sdJetCA.eta(),sdJetCA.phi(),sdConstituentsCA,1.5,20,20);

            tr.TriggerSubEvent("heat map");
          }

          //////////// subjet kinematics! /////////
          VJet *subjets = pfatjet->subjets;
          outjet->nsubjets=subjets->size();
          std::vector<sjpair> sjpairs;
          if (outjet->nsubjets>1) {
            // first set up the pairs
            double dR2 = DeltaR2(subjets->at(0),subjets->at(1));
            double mW = Mjj(subjets->at(0),subjets->at(1));
            double sumqg = subjets->at(0)->qgl + subjets->at(1)->qgl;
            sjpairs.emplace_back(dR2,mW,sumqg);
            outjet->sumqg=sumqg;
            outjet->minqg = TMath::Min(subjets->at(0)->qgl,subjets->at(1)->qgl);
            if (outjet->nsubjets>2) {
              dR2 = DeltaR2(subjets->at(0),subjets->at(2));
              mW = Mjj(subjets->at(0),subjets->at(2));
              sumqg = subjets->at(0)->qgl + subjets->at(2)->qgl;
              sjpairs.emplace_back(dR2,mW,sumqg);

              dR2 = DeltaR2(subjets->at(1),subjets->at(2));
              mW = Mjj(subjets->at(1),subjets->at(2));
              sumqg = subjets->at(1)->qgl + subjets->at(2)->qgl;
              sjpairs.emplace_back(dR2,mW,sumqg);

              // now order by dR
              std::sort(sjpairs.begin(),sjpairs.end(),orderByDR);
              outjet->dR2_minDR=sjpairs[0].dR2;
              outjet->mW_minDR=sjpairs[0].mW;

              // now by mW
              std::sort(sjpairs.begin(),sjpairs.end(),orderByMW);
              outjet->mW_best=sjpairs[0].mW;
              
              // now sumqg
              std::sort(sjpairs.begin(),sjpairs.end(),orderByQG);
              outjet->mW_qg=sjpairs[0].mW;
              outjet->sumqg += subjets->at(2)->qgl;
              outjet->minqg = TMath::Min(outjet->minqg,subjets->at(2)->qgl);
            } else {
              outjet->dR2_minDR=sjpairs[0].dR2;
              outjet->mW_minDR=subjets->at(0)->m;
              outjet->mW_best = (TMath::Abs(subjets->at(0)->m-WMASS)<TMath::Abs(subjets->at(1)->m-WMASS)) ? subjets->at(0)->m : subjets->at(1)->m;
              outjet->mW_qg=sjpairs[0].mW;
            }
            outjet->avgqg = outjet->sumqg/outjet->nsubjets;
          }

          tr.TriggerSubEvent("subjet kinematics");

          //////////// kinematic fit! ///////////
          if (doKinFit) {
            PJet *sj1=0, *sj2=0, *sjb=0;
            if (subjets->size()>=3) {
              VJet leadingSubjets(subjets->begin(),subjets->begin()+3);
              std::sort(leadingSubjets.begin(),leadingSubjets.end(),orderByCSV);
              PerformKinFit(fitter,fitresults,leadingSubjets[1],leadingSubjets[2],leadingSubjets[0]); 
              outjet->fitconv = (fitresults->converged) ? 1 : 0;
              if (fitresults->converged) {
                outjet->fitprob = fitresults->prob;
                outjet->fitchi2 = fitresults->chisq;
                outjet->fitmass = fitresults->fitmass;
                outjet->fitmassW = fitresults->fitmassW;
              }
            }
            tr.TriggerSubEvent("kinematic fit");

          }
          
          ////// pull angles! ///////
          unsigned int nS=sdsubjets.size();
          std::vector<TVector2> pulls; pulls.reserve(3);
          if (nS>0) {
            pulls.push_back(GetPull(sdsubjets[0]));
            outjet->betapull1 = GetPullBeta(pulls[0]);
            if (nS>1) {
              pulls.push_back(GetPull(sdsubjets[1]));
              outjet->betapull2 = GetPullBeta(pulls[1]);
              // 0:01, 1:02, 2:12
              outjet->alphapull1 = GetPullAlpha(sdsubjets[0],sdsubjets[1],pulls[0]);
              if (nS>2) {
                pulls.push_back(GetPull(sdsubjets[2]));
                outjet->betapull3 = GetPullBeta(pulls[2]);
                outjet->alphapull2 = GetPullAlpha(sdsubjets[0],sdsubjets[2],pulls[0]);
                outjet->alphapull3 = GetPullAlpha(sdsubjets[1],sdsubjets[2],pulls[1]);

                outjet->minpullangle = TMath::Abs(outjet->alphapull1); outjet->mW_minalphapull = Mjj(sdsubjets[0],sdsubjets[1]);
                if (TMath::Abs(outjet->alphapull2) < outjet->minpullangle) {
                  outjet->minpullangle = TMath::Abs(outjet->alphapull2);
                  outjet->mW_minalphapull = Mjj(sdsubjets[0],sdsubjets[2]);
                }
                if (TMath::Abs(outjet->alphapull3) < outjet->minpullangle) {
                  outjet->minpullangle = TMath::Abs(outjet->alphapull3);
                  outjet->mW_minalphapull = Mjj(sdsubjets[1],sdsubjets[2]);
                }
              }
              outjet->mW_minalphapull = sdsubjets[0].m();
            }
          }
          tr.TriggerSubEvent("pulls");

          /////// shower deco ///////
          // TODO 

          /////// fill ////////
          anafatjet->outtree->Fill();
        } else {
          //???
          PError("PandaAnalyzer::Run","No jet was clustered???");
        }

        tr.TriggerEvent("fat jet",false);
      } // loop over jets
    } // loop over jet collections


  } // entry loop

  if (DEBUG) { PDebug("PandaAnalyzer::Run","Done with entry loop"); }

} // Run()

