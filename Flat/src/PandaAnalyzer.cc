#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define DEBUG 0
using namespace panda;
using namespace std;

PandaAnalyzer::PandaAnalyzer() {
  gt = new GeneralTree();
  betas = gt->get_betas(); 
  Ns = gt->get_Ns(); 
  orders = gt->get_orders(); 
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
  TString jetname = (usePuppi) ? "puppi" : "chs";
  if (doFatjet)
    t->SetBranchAddress(jetname+"CA15",&fatjets);
  t->SetBranchAddress(jetname+"AK4",&jets);
  t->SetBranchAddress("electron",&electrons);
  t->SetBranchAddress("muon",&muons);
  t->SetBranchAddress("tau",&taus);
  t->SetBranchAddress("photon",&photons);
  t->SetBranchAddress("pfmet",&pfmet);
  t->SetBranchAddress("puppimet",&puppimet);
  if (!isData) {
    t->SetBranchAddress("gen",&genparts);  
    TH1F *hDTotalMCWeight = new TH1F("hDTotalMCWeight","hDTotalMCWeight",4,-2,2);
    TString val("fabs(info.mcWeight)/(info.mcWeight)");
    infotree->Draw(val+">>hDTotalMCWeight",val);
    fOut->WriteTObject(hDTotalMCWeight,"hDTotalMCWeight");
  }

}

PGenParticle *PandaAnalyzer::MatchToGen(double eta, double phi, double radius, int pdgid) {
  PGenParticle *found=NULL;
  double r2 = radius*radius;
  pdgid = abs(pdgid);

  unsigned int counter=0;
  for (map<PGenParticle*,float>::iterator iG=genObjects.begin();
        iG!=genObjects.end(); ++iG) {
    if (found!=NULL)
      break;
    if (pdgid!=0 && abs(iG->first->pdgid)!=pdgid)
      continue;
    if (DeltaR2(eta,phi,iG->first->eta,iG->first->phi)<r2) 
      found = iG->first;
  }

  return found;
}

void PandaAnalyzer::Terminate() {
  fOut->WriteTObject(tOut);
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
  if (preselBits & kMonotop) {
    if (gt->nFatjet>=1 && gt->fj1Pt>250) { 
      if ( (gt->puppimet>200 || gt->UZmag>200 || gt->UWmag>200 || gt->UAmag>200) ||
            (gt->pfmet>200 || gt->pfUZmag>200 || gt->pfUWmag>200 || gt->pfUAmag>200) ) {
            isGood = true;
      }
    }
  }
  if (preselBits & kMonojet) {
    if (gt->nJet>=1 && gt->jet1Pt>40) { 
      if ( (gt->puppimet>180 || gt->UZmag>180 || gt->UWmag>180 || gt->UAmag>180) ||
            (gt->pfmet>180 || gt->pfUZmag>180 || gt->pfUWmag>180 || gt->pfUAmag>180) ) {
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
      gt->sf_pu = getVal(hPUWeight,bound(gt->npv,sf_puMin,sf_puMax));
    if (isData) {
      for (auto iT : metTriggers) {   
        if (event->tiggers->at(iT)) {
          gt->trigger |= kMETTrig;
          break;
        }
      }
      for (auto iT : eleTriggers) {   
        if (event->tiggers->at(iT)) {
          gt->trigger |= kSingleEleTrig;
          break;
        }
      }
      for (auto iT : phoTriggers) {   
        if (event->tiggers->at(iT)) {
          gt->trigger |= kSinglePhoTrig;
          break;
        }
      }
    }


    tr.TriggerEvent("initialize");

    // met
    gt->pfmet = pfmet->pt;
    gt->pfmetphi = pfmet->phi;
    gt->puppimet = puppimet->pt;
    gt->puppimetphi = puppimet->phi;
    TLorentzVector vPFMET, vPuppiMET;
    vPFMET.SetPtEtaPhiM(gt->pfmet,0,gt->pfmetphi,0); 
    vPuppiMET.SetPtEtaPhiM(gt->puppimet,0,gt->puppimetphi,0); 
    TVector2 vMETNoMu; vMETNoMu.SetMagPhi(gt->pfmet,gt->pfmetphi); //  for trigger eff

    tr.TriggerEvent("met");

    //electrons
    std::vector<PObject*> looseLeps, tightLeps;
    for (PElectron *ele : *electrons) {
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
    for (PMuon *mu : *muons) {
      float pt = mu->pt; float eta = mu->eta; float aeta = fabs(eta);
      if (pt<10 || aeta>2.4)
        continue;
      if ((mu->id&PMuon::kLoose)==0)
        continue;
      if (!MuonIsolation(pt,eta,mu->iso,PMuon::kLoose))
        continue;
      looseLeps.push_back(mu);
      gt->nLooseMuon++;
      TVector2 vMu; vMu.SetMagPhi(pt,mu->phi);
      vMETNoMu += vMu;
    } 
    gt->pfmetnomu = vMETNoMu.Mod();

    // now consider all leptons
    gt->nLooseLep = looseLeps.size();
    if (gt->nLooseLep>0) {
      int nToSort = TMath::Min(3,gt->nLooseLep);
      std::partial_sort(looseLeps.begin(),looseLeps.begin()+nToSort,looseLeps.end(),SortPObjects);
    }
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
      PMuon *mu = dynamic_cast<PMuon*>(lep);
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
        PElectron *ele = dynamic_cast<PElectron*>(lep);
        bool isTight = ( (ele->id&PElectron::kTight)!=0 &&
                          ElectronIsolation(ele->pt,ele->eta,ele->iso,PElectron::kTight) &&
                          ele->pt>40 && fabs(ele->eta)<2.5 );
        if (lep_counter==1) {
          gt->looseLep1PdgId = ele->q*11;
          if (isTight) {
            gt->nTightElectron++;
            gt->looseLep1IsTight = 1;
            matchLeps.push_back(lep);
            matchEles.push_back(lep);
          }
        } else if (lep_counter==2) {
          gt->looseLep2PdgId = ele->q*11;
          if (isTight) {
            gt->nTightElectron++;
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
    if (gt->nLooseLep>1 && gt->looseLep1PdgId+gt->looseLep2PdgId==0) {
      TLorentzVector v1,v2;
      PObject *lep1=looseLeps[0], *lep2=looseLeps[1];
      v1.SetPtEtaPhiM(lep1->pt,lep1->eta,lep1->phi,lep1->m);
      v2.SetPtEtaPhiM(lep2->pt,lep2->eta,lep2->phi,lep2->m);
      gt->diLepMass = (v1+v2).M();
    } else {
      gt->diLepMass = -1;
    }

    tr.TriggerEvent("leptons"); 

    // photons
    std::vector<panda::PPhoton*> loosePhos;
    for (PPhoton *pho : *photons) {
      if ((pho->id&PPhoton::kLoose)==0)
        continue;
      float pt = pho->pt;
      if (pt<1) continue;
      float eta = pho->eta, phi = pho->phi;
      if (pt<15 || fabs(eta)>2.5)
        continue;
      if (IsMatched(&matchEles,0.16,eta,phi))
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

    // trigger efficiencies
    gt->sf_eleTrig=1; gt->sf_metTrig=1; gt->sf_phoTrig=1;
    if (!isData) {
      gt->sf_metTrig = getVal(hMetTrig,bound(gt->pfmetnomu,0,1000));

      if (gt->nLooseElectron>0 && abs(gt->looseLep1PdgId)==1
          && gt->looseLep1IsTight==1 && gt->nLooseMuon==0) {
        float eff1=0, eff2=0;
        if (gt->looseLep1Pt<100) {
          eff1 = getVal(hEleTrigLow,
                        bound(gt->looseLep1Eta,-1.*sf_eleEtaMax,sf_eleEtaMax),
                        bound(gt->looseLep1Pt,0,sf_elePtMax));
        } else {
          if (fabs(gt->looseLep1Eta)<1.4442) {
            eff1 = getVal(hEleTrigB,bound(gt->looseLep1Pt,0,1000));
          }
          if (1.5660<fabs(gt->looseLep1Eta) && fabs(gt->looseLep1Eta)<2.5) {
            eff1 = getVal(hEleTrigE,bound(gt->looseLep1Pt,0,1000));
          }
        }
        if (gt->nLooseElectron>1 && fabs(gt->looseLep2PdgId)==11) {
          if (gt->looseLep2Pt<100) {
            eff2 = getVal(hEleTrigLow,
                          bound(gt->looseLep2Eta,-1.*sf_eleEtaMax,sf_eleEtaMax),
                          bound(gt->looseLep2Pt,0,sf_elePtMax));
          } else {
            if (fabs(gt->looseLep2Eta)<1.4442) {
              eff2 = getVal(hEleTrigB,bound(gt->looseLep2Pt,0,1000));
            }
            if (1.5660<fabs(gt->looseLep2Eta) && fabs(gt->looseLep2Eta)<2.5) {
              eff2 = getVal(hEleTrigE,bound(gt->looseLep2Pt,0,1000));
            }
          }
        }
        gt->sf_eleTrig = 1 - (1-eff1)*(1-eff2);
      } // done with ele trig SF

      if (gt->nLoosePhoton>0 && gt->loosePho1IsTight)
        gt->sf_phoTrig = getVal(hPhoTrig,bound(gt->loosePho1Pt,160,1000));
    } 

    tr.TriggerEvent("triggers");

    // recoil!
    TLorentzVector vObj1, vObj2;
    TLorentzVector vUW, vUZ, vUA;
    TLorentzVector vpfUW, vpfUZ, vpfUA;
    if (gt->nLooseLep>0) {
      PObject *lep1 = looseLeps.at(0);
      vObj1.SetPtEtaPhiM(lep1->pt,lep1->eta,lep1->phi,lep1->m);

      // one lep => W
      vUW = vPuppiMET+vObj1; gt->UWmag=vUW.Pt(); gt->UWphi=vUW.Phi();
      vpfUW = vPFMET+vObj1; gt->pfUWmag=vpfUW.Pt(); gt->pfUWphi=vpfUW.Phi();

      if (gt->nLooseLep>1 && gt->looseLep1PdgId+gt->looseLep2PdgId==0) {
        // two OS lep => Z
        PObject *lep2 = looseLeps.at(1);
        vObj2.SetPtEtaPhiM(lep2->pt,lep2->eta,lep2->phi,lep2->m);

        vUZ=vUW+vObj2; gt->UZmag=vUZ.Pt(); gt->UZphi=vUZ.Phi();
        vpfUZ=vpfUW+vObj2; gt->pfUZmag=vpfUZ.Pt(); gt->pfUZphi=vpfUZ.Phi();
      }
    }
    if (gt->nLoosePhoton>0) {
      PPhoton *pho = loosePhos.at(0);
      vObj1.SetPtEtaPhiM(pho->pt,pho->eta,pho->phi,pho->m);

      vUA=vPuppiMET+vObj1; gt->UAmag=vUA.Pt(); gt->UAphi=vUA.Phi();
      vpfUA=vPFMET+vObj1; gt->pfUAmag=vpfUA.Pt(); gt->pfUAphi=vpfUA.Phi();
    }

    tr.TriggerEvent("recoils");

    PFatJet *fj1=0;
    gt->nFatjet=0;
    if (doFatjet) {
      for (PFatJet *fj : *fatjets) {
        float pt = fj->pt;
        float rawpt = fj->rawPt;
        float eta = fj->eta;
        float ptcut = 250;
        PDebug("",TString::Format("Considering fatjet with pt=%.3f, eta=%.3f",pt,eta));
        if (pt<ptcut || fabs(eta)>2.4)
          continue;
        float phi = fj->phi;
        if (IsMatched(&matchLeps,2.25,eta,phi) || IsMatched(&matchPhos,2.25,eta,phi)) {
          PDebug("",TString::Format("Rejecting fatjet with leps=%i, phos=%i",int(IsMatched(&matchLeps,2.25,eta,phi)),int(IsMatched(&matchPhos,2.25,eta,phi))));
          continue;
          /*
          if (gt->nFatjet==0) {
            break;
          } else {
            continue;
          }
          */
        }
        gt->nFatjet++;
        if (gt->nFatjet==1) {
          fj1 = fj;
          gt->fj1Pt = pt;
          gt->fj1Eta = eta;
          gt->fj1Phi = phi;
          gt->fj1MSD = fj->mSD;
          gt->fj1Tau32 = clean(fj->tau3/fj->tau2);
          gt->fj1Tau32SD = clean(fj->tau3SD/fj->tau2SD);
          gt->fj1Tau21 = clean(fj->tau2/fj->tau1);
          gt->fj1Tau21SD = clean(fj->tau2SD/fj->tau1SD);
          gt->fj1RawPt = rawpt;

          for (unsigned int iB=0; iB!=betas.size(); ++iB) {
            float beta = betas.at(iB);
            for (auto N : Ns) {
              for (auto order : orders) {
                gt->fj1ECFNs[makeECFString(order,N,beta)] = fj->get_ecf(order,N,iB); 
              }
            }
          } //loop over betas

          VJet *subjets = fj->subjets;
          std::sort(subjets->begin(),subjets->end(),SortPJetByCSV);
          gt->fj1MaxCSV = subjets->at(0)->csv; 
          gt->fj1MinCSV = subjets->back()->csv; 
        }
      }
    }

    tr.TriggerEvent("fatjet");

    // first identify interesting jets
    vector<PJet*> cleanedJets, isoJets;
    TLorentzVector vJet;
    PJet *jet1=0, *jet2=0;
    gt->dphipuppimet=999; gt->dphipfmet=999;
    for (PJet *jet : *jets) {
      if (jet->pt<30 || abs(jet->eta)>4.5) // loose ID should go here
        continue;
      if (IsMatched(&matchLeps,0.16,jet->eta,jet->phi) ||
          IsMatched(&matchPhos,0.16,jet->eta,jet->phi))
        continue;
      cleanedJets.push_back(jet);
      float csv = (fabs(jet->eta)<2.5) ? jet->csv : -1;
      if (cleanedJets.size()==1) {
        jet1 = jet;
        gt->jet1Pt = jet->pt;
        gt->jet1Eta = jet->eta;
        gt->jet1Phi = jet->phi;
        gt->jet1CSV = csv; 
      } else if (cleanedJets.size()==2) {
        jet2 = jet;
        gt->jet2Pt = jet->pt;
        gt->jet2Eta = jet->eta;
        gt->jet2Phi = jet->phi;
        gt->jet2CSV = csv; 
      } 
      // compute dphi wrt mets
      vJet.SetPtEtaPhiM(jet->pt,jet->eta,jet->phi,jet->m);
      gt->dphipuppimet = std::min(fabs(vJet.DeltaPhi(vPuppiMET)),(double)gt->dphipuppimet);
      gt->dphipfmet = std::min(fabs(vJet.DeltaPhi(vPFMET)),(double)gt->dphipfmet);
      // btags
      if (csv>0.460) ++gt->jetNBtags;
      if (gt->nFatjet>0 && fabs(jet->eta)<2.5
          && DeltaR2(gt->fj1Eta,gt->fj1Phi,jet->eta,jet->phi)>2.25) {
        isoJets.push_back(jet);
        if (csv>0.460) ++gt->isojetNBtags;
      }
    }
    gt->nJet = cleanedJets.size();
    if (gt->nJet>1) {
      gt->jet12DEta = fabs(jet1->eta-jet2->eta);
      TLorentzVector vj1, vj2; 
      vj1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->m);
      vj2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->m);
      gt->jet12Mass = (vj1+vj2).M();
    }

    tr.TriggerEvent("jets");

    for (PTau *tau : *taus) {
      if ((tau->id&PTau::kDecayModeFinding)==0 ||
          (tau->id&PTau::kDecayModeFindingNewDMs)==0) 
        continue;
      if (tau->isoDeltaBetaCorr>5)
        continue; 
      if (tau->pt<10 || fabs(tau->eta)>2.3)
        continue;
      if (IsMatched(&matchLeps,0.16,tau->eta,tau->phi))
        continue;
      gt->nTau++;
    } 

    tr.TriggerEvent("taus");

    if (!PassPreselection())
      continue;

    tr.TriggerEvent("presel");

    // identify interesting gen particles for fatjet matching
    unsigned int pdgidTarget=0;
    if (!isData && processType>=kTT) {
      switch(processType) {
        case kTop:
        case kTT:
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

      int nGen = genparts->size();
      for (int iG=0; iG!=nGen; ++iG) {
        PGenParticle *part = genparts->at(iG);
        int pdgid = part->pdgid;
        unsigned int abspdgid = abs(pdgid);
        if (abspdgid == pdgidTarget)
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
        if (processType==kTop||processType==kTT) {
          
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
          if (iW<0) {// ???
            PDebug("","Could not find W");
            continue;
          }
          PGenParticle *partW = genparts->at(iW);

          // now look for b or W->qq
          int iB=-1, iQ1=-1, iQ2=-1;
          double size=0, sizeW=0;
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
                sizeW = TMath::Max(DeltaR2(partW->eta,partW->phi,partQ->eta,partQ->phi),sizeW);
              } else if (iQ2<0) {
                iQ2 = jG;
                size = TMath::Max(DeltaR2(part->eta,part->phi,partQ->eta,partQ->phi),size);
                sizeW = TMath::Max(DeltaR2(partW->eta,partW->phi,partQ->eta,partQ->phi),sizeW);
              }
            }
            if (iB>=0 && iQ1>=0 && iQ2>=0)
              break;
          } // looking for quarks


          bool isHadronic = (iB>=0 && iQ1>=0 && iQ2>=0); // all 3 quarks were found
          if (isHadronic)
            genObjects[part] = size;

          bool isHadronicW = (iQ1>=0 && iQ2>=0);
          if (isHadronicW)
            genObjects[partW] = sizeW;

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

    // do gen matching now that presel is passed
    if (!isData && gt->nFatjet>0) {
      // first see if jet is matched
      PGenParticle *matched = MatchToGen(fj1->eta,fj1->phi,1.5,pdgidTarget);
      if (matched!=NULL) {
        gt->fj1IsMatched = 1;
        gt->fj1GenPt = matched->pt;
        gt->fj1GenSize = genObjects[matched];
      } else {
        gt->fj1IsMatched = 0;
      }
      if (pdgidTarget==6) { // matched to top; try for W
        PGenParticle *matchedW = MatchToGen(fj1->eta,fj1->phi,1.5,24);
        if (matchedW!=NULL) {
          gt->fj1IsWMatched = 1;
          gt->fj1GenWPt = matchedW->pt;
          gt->fj1GenWSize = genObjects[matchedW];
        } else {
          gt->fj1IsWMatched = 0;
        }
      }
      
      // now get the highest pT gen particle inside the jet cone
      for (PGenParticle *gen : *genparts) {
        float pt = gen->pt;
        int pdgid = gen->pdgid;
        if (pt>(gt->fj1HighestPtGenPt) 
            && DeltaR2(gen->eta,gen->phi,fj1->eta,fj1->phi)<2.25) {
          gt->fj1HighestPtGenPt = pt;
          gt->fj1HighestPtGen = pdgid;
        }
      }

      // now get the subjet btag SFs
      vector<btagcand> sj_btagcands;
      vector<double> sj_sf_cent, sj_sf_bUp, sj_sf_bDown, sj_sf_mUp, sj_sf_mDown;
      unsigned int nSJ = fj1->subjets->size();
      for (unsigned int iSJ=0; iSJ!=nSJ; ++iSJ) {
        PJet *subjet = fj1->subjets->at(iSJ);
        int flavor=0;
        for (PGenParticle *gen : *genparts) {
          int apdgid = abs(gen->pdgid);
          if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
            continue;
          double dr2 = DeltaR2(subjet->eta,subjet->phi,gen->eta,gen->phi);
          if (dr2<0.09) {
            if (apdgid==4 || apdgid==5) {
              flavor=apdgid;
              break;
            } else {
              flavor=0;
            }
          }
        } // finding the subjet flavor
        float sjPtMax = (flavor<4) ? 1000. : 450.;
        float pt = subjet->pt;
        float btagUncFactor = 1;
        if (pt>sjPtMax) {
          btagUncFactor = 2;
          pt = sjPtMax;
        }
        float eta = subjet->eta;
        double eff,sf,sfUp,sfDown;
        unsigned int bin = btagpt.bin(pt);
        if (flavor==5) {
          eff = beff[bin];
          sf = sj_hfReader->eval(BTagEntry::FLAV_B,eta,pt,0);
          sfUp = sj_hfUpReader->eval(BTagEntry::FLAV_B,eta,pt,0);
          sfDown = sj_hfDownReader->eval(BTagEntry::FLAV_B,eta,pt,0);
        } else if (flavor==4) {
          eff = ceff[bin];
          sf = sj_hfReader->eval(BTagEntry::FLAV_C,eta,pt,0);
          sfUp = sj_hfUpReader->eval(BTagEntry::FLAV_C,eta,pt,0);
          sfDown = sj_hfDownReader->eval(BTagEntry::FLAV_C,eta,pt,0);
        } else {
          eff = lfeff[bin];
          sf = sj_lfReader->eval(BTagEntry::FLAV_UDSG,eta,pt,0);
          sfUp = sj_lfUpReader->eval(BTagEntry::FLAV_UDSG,eta,pt,0);
          sfDown = sj_lfDownReader->eval(BTagEntry::FLAV_UDSG,eta,pt,0);
        }
        sfUp = btagUncFactor*(sfUp-sf)+sf;
        sfDown = btagUncFactor*(sfDown-sf)+sf;
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
      EvalBtagSF(sj_btagcands,sj_sf_cent,
                  gt->sf_sjbtag0,gt->sf_sjbtag1,gt->sf_sjbtag2);
      EvalBtagSF(sj_btagcands,sj_sf_bUp,
                  gt->sf_sjbtag0BUp,gt->sf_sjbtag1BUp,gt->sf_sjbtag2BUp);
      EvalBtagSF(sj_btagcands,sj_sf_bDown,
                  gt->sf_sjbtag0BDown,gt->sf_sjbtag1BDown,gt->sf_sjbtag2BDown);
      EvalBtagSF(sj_btagcands,sj_sf_mUp,
                  gt->sf_sjbtag0MUp,gt->sf_sjbtag1MUp,gt->sf_sjbtag2MUp);
      EvalBtagSF(sj_btagcands,sj_sf_mDown,
                  gt->sf_sjbtag0MDown,gt->sf_sjbtag1MDown,gt->sf_sjbtag2MDown);

    }

    tr.TriggerEvent("fatjet gen-matching");

    if (!isData) {
      // now get the jet btag SFs
      vector<btagcand> btagcands;
      vector<double> sf_cent, sf_bUp, sf_bDown, sf_mUp, sf_mDown;
      unsigned int nJ = isoJets.size();
      for (unsigned int iJ=0; iJ!=nJ; ++iJ) {
        PJet *jet = isoJets.at(iJ);
        int flavor=0;
        for (PGenParticle *gen : *genparts) {
          int apdgid = abs(gen->pdgid);
          if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
            continue;
          double dr2 = DeltaR2(jet->eta,jet->phi,gen->eta,gen->phi);
          if (dr2<0.09) {
            if (apdgid==4 || apdgid==5) {
              flavor=apdgid;
              break;
            } else {
              flavor=0;
            }
          }
        } // finding the jet flavor
        float jPtMax = 670.;
        float pt = jet->pt;
        float btagUncFactor = 1;
        if (pt>jPtMax) {
          btagUncFactor = 2;
          pt = jPtMax;
        }
        float eta = jet->eta;
        double eff,sf,sfUp,sfDown;
        unsigned int bin = btagpt.bin(pt);
        if (flavor==5) {
          eff = beff[bin];
          sf = hfReader->eval(BTagEntry::FLAV_B,eta,pt,0);
          sfUp = hfUpReader->eval(BTagEntry::FLAV_B,eta,pt,0);
          sfDown = hfDownReader->eval(BTagEntry::FLAV_B,eta,pt,0);
        } else if (flavor==4) {
          eff = ceff[bin];
          sf = hfReader->eval(BTagEntry::FLAV_C,eta,pt,0);
          sfUp = hfUpReader->eval(BTagEntry::FLAV_C,eta,pt,0);
          sfDown = hfDownReader->eval(BTagEntry::FLAV_C,eta,pt,0);
        } else {
          eff = lfeff[bin];
          sf = lfReader->eval(BTagEntry::FLAV_UDSG,eta,pt,0);
          sfUp = lfUpReader->eval(BTagEntry::FLAV_UDSG,eta,pt,0);
          sfDown = lfDownReader->eval(BTagEntry::FLAV_UDSG,eta,pt,0);
        }
        sfUp = btagUncFactor*(sfUp-sf)+sf;
        sfDown = btagUncFactor*(sfDown-sf)+sf;
        btagcands.push_back(btagcand(iJ,flavor,eff,sf,sfUp,sfDown));
        sf_cent.push_back(sf);
        if (flavor>0) {
          sf_bUp.push_back(sfUp); sf_bDown.push_back(sfDown);
          sf_mUp.push_back(sf); sf_mDown.push_back(sf);
        } else {
          sf_bUp.push_back(sf); sf_bDown.push_back(sf);
          sf_mUp.push_back(sfUp); sf_mDown.push_back(sfDown);
        }
      } // loop over jets
      EvalBtagSF(btagcands,sf_cent,
                  gt->sf_btag0,gt->sf_btag1);
      EvalBtagSF(btagcands,sf_bUp,
                  gt->sf_btag0BUp,gt->sf_btag1BUp);
      EvalBtagSF(btagcands,sf_bDown,
                  gt->sf_btag0BDown,gt->sf_btag1BDown);
      EvalBtagSF(btagcands,sf_mUp,
                  gt->sf_btag0MUp,gt->sf_btag1MUp);
      EvalBtagSF(btagcands,sf_mDown,
                  gt->sf_btag0MDown,gt->sf_btag1MDown);
    }

    tr.TriggerEvent("ak4 gen-matching");

    // ttbar pT weight
    gt->sf_tt = 1;
    gt->sf_tt_ext = 1;
    if (!isData && processType==kTT) {
      float pt_t=0, pt_tbar=0;
      for (auto *gen : *genparts) {
        if (abs(gen->pdgid)!=6)
          continue;
        if (gen->parent>=0 && genparts->at(gen->parent)->pdgid==gen->pdgid)
          continue; // must be first copy
        if (gen->pdgid>0)
          pt_t = gen->pt;
        else
          pt_tbar = gen->pt;
        if (pt_t>0 && pt_tbar>0)
          break;
      }
      if (pt_t>0 && pt_tbar>0) {
        gt->sf_tt_ext = TMath::Sqrt( TMath::Exp(0.156-0.00137*pt_t) * TMath::Exp(0.156-0.00137*pt_tbar) ); // extend the fit past 400
        pt_t = bound(pt_t,0,400);
        pt_tbar = bound(pt_tbar,0,400);
        gt->sf_tt = TMath::Sqrt( TMath::Exp(0.156-0.00137*pt_t) * TMath::Exp(0.156-0.00137*pt_tbar) ); // bound at 400
      }
    } 

    tr.TriggerEvent("tt SFs");

    // derive ewk/qcd weights
    gt->sf_qcdV=1; gt->sf_ewkV=1;
    if (!isData) {
      bool found = processType!=kA && processType!=kZ && processType!=kW;
      int target=24;
      if (processType==kZ) target=23;
      if (processType==kA) target=22;
      for (PGenParticle *gen : *genparts) {
        if (found) break;
        int apdgid = abs(gen->pdgid);
        if (apdgid==target)  {
          if (gen->parent>=0 && genparts->at(gen->parent)->pdgid==gen->pdgid)
            continue;
          if (processType==kZ) {
            gt->trueGenBosonPt = gen->pt;
            gt->genBosonPt = bound(gen->pt,genBosonPtMin,genBosonPtMax);
            gt->sf_qcdV = getVal(hZNLO,gt->genBosonPt);
            gt->sf_ewkV = getVal(hZEWK,gt->genBosonPt);
            found=true;
          } else if (processType==kW) {
            gt->trueGenBosonPt = gen->pt;
            gt->genBosonPt = bound(gen->pt,genBosonPtMin,genBosonPtMax);
            gt->sf_qcdV = getVal(hWNLO,gt->genBosonPt);
            gt->sf_ewkV = getVal(hWEWK,gt->genBosonPt);
            found=true;
          } else if (processType==kA) {
            // take the highest pT
            if (gen->pt > gt->trueGenBosonPt) {
              gt->trueGenBosonPt = gen->pt;
              gt->genBosonPt = bound(gen->pt,genBosonPtMin,genBosonPtMax);
              gt->sf_qcdV = getVal(hANLO,gt->genBosonPt);
              gt->sf_ewkV = getVal(hAEWK,gt->genBosonPt);
            }
          }
        } // target matches
      }
    }

    tr.TriggerEvent("qcd/ewk SFs");

    //lepton SFs
    gt->sf_lep=1; gt->sf_lepTrack=1;
    if (!isData) {
      for (unsigned int iL=0; iL!=TMath::Min(gt->nLooseLep,2); ++iL) {
        PObject *lep = looseLeps.at(iL);
        float pt = lep->pt, eta = lep->eta;
        bool isTight = (iL==0 && gt->looseLep1IsTight) || (iL==1 && gt->looseLep2IsTight);
        PMuon *mu = dynamic_cast<PMuon*>(lep);
        if (mu!=NULL) {
          pt = bound(pt,sf_muPtMin,sf_muPtMax);
          eta = bound(eta,0,sf_muEtaMax);
          if (isTight)
            gt->sf_lep *= getVal(hMuTight,eta,pt);
          else
            gt->sf_lep *= getVal(hMuLoose,eta,pt);
          gt->sf_lepTrack *= getVal(hMuTrack,bound(gt->npv,0,50));
        } else {
          PElectron *ele = dynamic_cast<PElectron*>(lep);
          pt = bound(pt,sf_elePtMin,sf_elePtMax);
          eta = bound(eta,0,sf_eleEtaMax);
          if (isTight)
            gt->sf_lep *= getVal(hEleTight,eta,pt);
          else
            gt->sf_lep *= getVal(hEleVeto,eta,pt);
          gt->sf_lepTrack *= getVal(hEleTrack,bound(ele->eta,-2.5,2.5),bound(gt->npv,0,50));
        }
      }
    }

    tr.TriggerEvent("lepton SFs");

    gt->Fill();

  } // entry loop

  if (DEBUG) { PDebug("PandaAnalyzer::Run","Done with entry loop"); }

} // Run()

