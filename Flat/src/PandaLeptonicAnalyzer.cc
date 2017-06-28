#include "../interface/PandaLeptonicAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#include "TSystem.h"
#include "TRandom2.h"

#include "PandaAnalysis/Flat/interface/LeptonicFactors.h"

using namespace panda;
using namespace std;

PandaLeptonicAnalyzer::PandaLeptonicAnalyzer(int debug_/*=0*/) {
  DEBUG = debug_;

  if (DEBUG) PDebug("PandaLeptonicAnalyzer::PandaLeptonicAnalyzer","Calling constructor");
  gt = new GeneralLeptonicTree();
  if (DEBUG) PDebug("PandaLeptonicAnalyzer::PandaLeptonicAnalyzer","Built GeneralLeptonicTree");
  flags["puppi"]     = false;
  flags["applyJSON"] = true;
  flags["genOnly"]   = false;
  flags["lepton"]    = false;
  if (DEBUG) PDebug("PandaLeptonicAnalyzer::PandaLeptonicAnalyzer","Called constructor");
}


PandaLeptonicAnalyzer::~PandaLeptonicAnalyzer() {
  if (DEBUG) PDebug("PandaLeptonicAnalyzer::~PandaLeptonicAnalyzer","Calling destructor");
}


void PandaLeptonicAnalyzer::ResetBranches() {
  genObjects.clear();
  matchPhos.clear();
  matchEles.clear();
  matchLeps.clear();
  gt->Reset();
  if (DEBUG) PDebug("PandaLeptonicAnalyzer::ResetBranches","Reset");
}


void PandaLeptonicAnalyzer::SetOutputFile(TString fOutName) {
  fOut = new TFile(fOutName,"RECREATE");
  tOut = new TTree("events","events");

  fOut->WriteTObject(hDTotalMCWeight);

  // fill the signal weights
  for (auto& id : wIDs) 
    gt->signal_weights[id] = 1;

  // Build the input tree here 
  gt->WriteTree(tOut);

  if (DEBUG) PDebug("PandaLeptonicAnalyzer::SetOutputFile","Created output in "+fOutName);
}


int PandaLeptonicAnalyzer::Init(TTree *t, TH1D *hweights, TTree *weightNames)
{
  if (DEBUG) PDebug("PandaLeptonicAnalyzer::Init","Starting initialization");
  if (!t || !hweights) {
    PError("PandaLeptonicAnalyzer::Init","Malformed input!");
    return 0;
  }
  tIn = t;

  event.setStatus(*t, {"!*"}); // turn everything off first

  TString jetname = (flags["puppi"]) ? "puppi" : "chs";
  panda::utils::BranchList readlist({"runNumber", "lumiNumber", "eventNumber", "rho", 
                                     "isData", "npv", "npvTrue", "weight", "chsAK4Jets", 
                                     "electrons", "muons", "taus", "photons", 
                                     "pfMet", "caloMet", "puppiMet", "rawMet", 
                                     "recoil","metFilters","genMet",});
  readlist.setVerbosity(0);

  readlist.push_back("triggers");

  if (isData) {
  } else {
   readlist.push_back("genParticles");
   readlist.push_back("genReweight");
  }

  event.setAddress(*t, readlist); // pass the readlist so only the relevant branches are turned on
  if (DEBUG) PDebug("PandaLeptonicAnalyzer::Init","Set addresses");

  hDTotalMCWeight = new TH1F("hDTotalMCWeight","hDTotalMCWeight",1,0,2);
  hDTotalMCWeight->SetBinContent(1,hweights->GetBinContent(1));

  const int nBinPt = 57; Float_t xbinsPt[nBinPt+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                       10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                       20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                                                       30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                                                       40, 50, 60, 70, 80, 90,100,125,150,175, 
						      200,250,300,350,400,450,500,1000};
  hDDilPtMM = new TH1D("hDDilPtMM", "hDDilPtMM", nBinPt, xbinsPt);
  hDDilPtEE = new TH1D("hDDilPtEE", "hDDilPtEE", nBinPt, xbinsPt);

  hDDilLowPtMM = new TH1D("hDDilLowPtMM", "hDDilLowPtMM", 500, 0, 50);
  hDDilLowPtEE = new TH1D("hDDilLowPtEE", "hDDilLowPtEE", 500, 0, 50);

  hDDilPt2MM = new TH1D("hDDilPt2MM", "hDDilPt2MM", 50, 0, 50);
  hDDilPt2EE = new TH1D("hDDilPt2EE", "hDDilPt2EE", 50, 0, 50);

  if (weightNames) {
    if (weightNames->GetEntries()!=377 && weightNames->GetEntries()!=22) {
      PError("PandaLeptonicAnalyzer::Init",
          TString::Format("Reweighting failed because only found %u weights!",
                          unsigned(weightNames->GetEntries())));
      return 1;
    }
    TString *id = new TString();
    weightNames->SetBranchAddress("id",&id);
    unsigned nW = weightNames->GetEntriesFast();
    for (unsigned iW=0; iW!=nW; ++iW) {
      weightNames->GetEntry(iW);
      wIDs.push_back(*id);
    }
  } else if (processType==kSignal) {
    PError("PandaLeptonicAnalyzer::Init","This is a signal file, but the weights are missing!");
    return 2;
  }


  // manipulate the output tree
  if (isData) {
    std::vector<TString> droppable = {"mcWeight","scale","pdf.*","gen.*","sf_.*"};
    gt->RemoveBranches(droppable,{"sf_phoPurity"});
  }
  if (flags["lepton"]) {
    std::vector<TString> droppable = {"sf_sjbtag*"};
    gt->RemoveBranches(droppable);
  }
  if (flags["genOnly"]) {
    std::vector<TString> keepable = {"mcWeight","scale","scaleUp",
                                     "scaleDown","pdf*","gen*",
                                     "sf_tt*","sf_qcdTT*"};
    gt->RemoveBranches({".*"},keepable);
  }

  if (DEBUG) PDebug("PandaLeptonicAnalyzer::Init","Finished configuration");

  return 0;
}


panda::GenParticle const *PandaLeptonicAnalyzer::MatchToGen(double eta, double phi, double radius, int pdgid) {
  panda::GenParticle const* found=NULL;
  double r2 = radius*radius;
  pdgid = abs(pdgid);

  unsigned int counter=0;
  for (map<panda::GenParticle const*,float>::iterator iG=genObjects.begin();
      iG!=genObjects.end(); ++iG) {
    if (found!=NULL)
      break;
    if (pdgid!=0 && abs(iG->first->pdgid)!=pdgid)
      continue;
    if (DeltaR2(eta,phi,iG->first->eta(),iG->first->phi())<r2)
      found = iG->first;
  }

  return found;
}


void PandaLeptonicAnalyzer::Terminate() {
  fOut->WriteTObject(tOut);
  fOut->WriteTObject(hDDilPtMM);    
  fOut->WriteTObject(hDDilPtEE);    
  fOut->WriteTObject(hDDilLowPtMM);    
  fOut->WriteTObject(hDDilLowPtEE);    
  fOut->WriteTObject(hDDilPt2MM);    
  fOut->WriteTObject(hDDilPt2EE);    
  fOut->Close();

  //for (auto *f : fCorrs)
  //  if (f)
  //    f->Close();
  for (auto *h : h1Corrs)
    delete h;
  for (auto *h : h2Corrs)
    delete h;

  delete btagCalib;
  for (auto *reader : btagReaders )
    delete reader;

  for (auto& iter : ak8UncReader)
    delete iter.second;

  delete ak8JERReader;

  for (auto& iter : ak4UncReader)
    delete iter.second;

  for (auto& iter : ak4ScaleReader)
    delete iter.second;

  delete ak4JERReader;
  
  delete hDTotalMCWeight;
  delete hDDilPtMM;
  delete hDDilPtEE;
  delete hDDilLowPtMM;
  delete hDDilLowPtEE;
  delete hDDilPt2MM;
  delete hDDilPt2EE;

  if (DEBUG) PDebug("PandaLeptonicAnalyzer::Terminate","Finished with output");
}

void PandaLeptonicAnalyzer::OpenCorrection(CorrectionType ct, TString fpath, TString hname, int dim) {
  fCorrs[ct] = TFile::Open(fpath);
  if (dim==1) 
    h1Corrs[ct] = new THCorr1((TH1D*)fCorrs[ct]->Get(hname));
  else
    h2Corrs[ct] = new THCorr2((TH2D*)fCorrs[ct]->Get(hname));
}

double PandaLeptonicAnalyzer::GetCorr(CorrectionType ct, double x, double y) {
  if (h1Corrs[ct]!=0) {
    return h1Corrs[ct]->Eval(x); 
  } else if (h2Corrs[ct]!=0) {
    return h2Corrs[ct]->Eval(x,y);
  } else {
    PError("PandaLeptonicAnalyzer::GetCorr",
       TString::Format("No correction is defined for CorrectionType=%u",ct));
    return 1;
  }
}

void PandaLeptonicAnalyzer::SetDataDir(const char *s2) {
  TString dirPath1 = TString(gSystem->Getenv("CMSSW_BASE")) + "/src/";
  TString dirPath2(s2);

  if (DEBUG) PDebug("PandaLeptonicAnalyzer::SetDataDir","Starting loading of data");

  // pileup
  OpenCorrection(cPU    ,dirPath1+"MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root","puWeights",1);
  OpenCorrection(cPUUp  ,dirPath1+"MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root","puWeightsUp",1);
  OpenCorrection(cPUDown,dirPath1+"MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root","puWeightsDown",1);

  OpenCorrection(ZHEwkCorr    ,dirPath1+"MitAnalysisRunII/data/80x/Zll_nloEWK_weight_unnormalized.root","SignalWeight_nloEWK_rebin",1);
  OpenCorrection(ZHEwkCorrUp  ,dirPath1+"MitAnalysisRunII/data/80x/Zll_nloEWK_weight_unnormalized.root","SignalWeight_nloEWK_up_rebin",1);
  OpenCorrection(ZHEwkCorrDown,dirPath1+"MitAnalysisRunII/data/80x/Zll_nloEWK_weight_unnormalized.root","SignalWeight_nloEWK_down_rebin",1);

  OpenCorrection(cLooseMuonId  ,dirPath1+"MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root","scalefactors_MuonLooseId_Muon",2);
  OpenCorrection(cMediumMuonId ,dirPath1+"MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root","scalefactors_MuonMediumId_Muon",2);
  OpenCorrection(cTightMuonId  ,dirPath1+"MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root","scalefactors_TightId_Muon",2);
  OpenCorrection(cLooseMuonIso ,dirPath1+"MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root","scalefactors_Iso_MuonLooseId",2);
  OpenCorrection(cMediumMuonIso,dirPath1+"MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root","scalefactors_Iso_MuonMediumId",2);
  OpenCorrection(cTightMuonIso ,dirPath1+"MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root","scalefactors_Iso_MuonTightId",2);
  OpenCorrection(cTrackingMuon ,dirPath1+"MitAnalysisRunII/data/80x/Tracking_EfficienciesAndSF_BCDEFGH.root","ratio_eff_eta3_dr030e030_corr",1);

  OpenCorrection(cLooseElectronId ,dirPath1+"MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root","scalefactors_Loose_Electron",2);
  OpenCorrection(cMediumElectronId,dirPath1+"MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root","scalefactors_Medium_Electron",2);
  OpenCorrection(cTightElectronId ,dirPath1+"MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root","scalefactors_Tight_Electron",2);
  OpenCorrection(cTrackingElectron,dirPath1+"MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root","scalefactors_Reco_Electron",2);

  // btag SFs
  btagCalib = new BTagCalibration("csvv2",(dirPath1+"MitAnalysisRunII/data/80x/CSVv2_Moriond17_B_H.csv").Data());
  btagReaders[bJetL] = new BTagCalibrationReader(BTagEntry::OP_LOOSE,"central",{"up","down"});
  btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_B,"comb");
  btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_C,"comb");
  btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_UDSG,"incl");

  btagReaders[bJetM] = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"});
  btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_B,"comb");
  btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_C,"comb");
  btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_UDSG,"incl");

  if (DEBUG) PDebug("PandaLeptonicAnalyzer::SetDataDir","Loaded btag SFs");

  TString jecV = "V4", jecReco = "23Sep2016"; 
  TString jecVFull = jecReco+jecV;
  ak8UncReader["MC"] = new JetCorrectionUncertainty(
     (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_Uncertainty_AK8PFPuppi.txt").Data()
    );
  std::vector<TString> eraGroups = {"BCD","EF","G","H"};
  for (auto e : eraGroups) {
    ak8UncReader["data"+e] = new JetCorrectionUncertainty(
       (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_Uncertainty_AK8PFPuppi.txt").Data()
      );
  }

  ak8JERReader = new JERReader(dirPath2+"/jec/25nsV10/Spring16_25nsV10_MC_SF_AK8PFPuppi.txt",
                               dirPath2+"/jec/25nsV10/Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt");


  ak4UncReader["MC"] = new JetCorrectionUncertainty(
     (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_Uncertainty_AK4PFPuppi.txt").Data()
    );
  for (auto e : eraGroups) {
    ak4UncReader["data"+e] = new JetCorrectionUncertainty(
       (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_Uncertainty_AK4PFPuppi.txt").Data()
      );
  }

  ak4JERReader = new JERReader(dirPath2+"/jec/25nsV10/Spring16_25nsV10_MC_SF_AK4PFPuppi.txt",
                               dirPath2+"/jec/25nsV10/Spring16_25nsV10_MC_PtResolution_AK4PFPuppi.txt");

  std::vector<JetCorrectorParameters> params = {
    JetCorrectorParameters(
      (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L1FastJet_AK4PFPuppi.txt").Data()),
    JetCorrectorParameters(
      (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L2Relative_AK4PFPuppi.txt").Data()),
    JetCorrectorParameters(
      (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L3Absolute_AK4PFPuppi.txt").Data()),
    JetCorrectorParameters(
      (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L2L3Residual_AK4PFPuppi.txt").Data())
  };
  ak4ScaleReader["MC"] = new FactorizedJetCorrector(params);
  if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::SetDataDir","Loaded JES for AK4 MC");
  for (auto e : eraGroups) {
    params = {
      JetCorrectorParameters(
        (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L1FastJet_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L2Relative_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L3Absolute_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L2L3Residual_AK4PFPuppi.txt").Data())
    };
    ak4ScaleReader["data"+e] = new FactorizedJetCorrector(params);
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::SetDataDir","Loaded JES for AK4 "+e);
  }

  if (DEBUG) PDebug("PandaLeptonicAnalyzer::SetDataDir","Loaded JES/R");

}

void PandaLeptonicAnalyzer::AddGoodLumiRange(int run, int l0, int l1) {
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) { // don't know about this run yet
    std::vector<LumiRange> newLumiList;
    newLumiList.emplace_back(l0,l1);
    goodLumis[run] = newLumiList;
  } else {
    run_->second.emplace_back(l0,l1);
  }
}


bool PandaLeptonicAnalyzer::PassGoodLumis(int run, int lumi) {
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) {
    // matched no run
    if (DEBUG) 
      PDebug("PandaLeptonicAnalyzer::PassGoodLumis",TString::Format("Failing run=%i",run));
    return false;
  }

  // found the run, now look for a lumi range
  for (auto &range : run_->second) {
    if (range.Contains(lumi)) {
      if (DEBUG) 
        PDebug("PandaLeptonicAnalyzer::PassGoodLumis",TString::Format("Accepting run=%i, lumi=%i",run,lumi));
      return true;
    }
  }

  // matched no lumi range
  if (DEBUG) 
    PDebug("PandaLeptonicAnalyzer::PassGoodLumis",TString::Format("Failing run=%i, lumi=%i",run,lumi));
  return false;
}


bool PandaLeptonicAnalyzer::PassPreselection() {
  
  if (preselBits==0)
    return true;
  bool isGood=false;

  if (preselBits & kLepton) {
    if (gt->nLooseLep>1 && gt->looseLep1Pt > 20 && gt->looseLep2Pt > 20) isGood = true;
  }

  return isGood;
}


void PandaLeptonicAnalyzer::CalcBJetSFs(BTagType bt, int flavor,
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

void PandaLeptonicAnalyzer::EvalBTagSF(std::vector<btagcand> &cands, std::vector<double> &sfs,
               GeneralLeptonicTree::BTagShift shift,GeneralLeptonicTree::BTagJet jettype, bool do2) 
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

  GeneralLeptonicTree::BTagParams p;
  p.shift = shift;
  p.jet = jettype;
  p.tag=GeneralLeptonicTree::b0; gt->sf_btags[p] = sf0;
  p.tag=GeneralLeptonicTree::b1; gt->sf_btags[p] = sf1;
  p.tag=GeneralLeptonicTree::bGT0; gt->sf_btags[p] = sfGT0;

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

    p.tag=GeneralLeptonicTree::b2; gt->sf_btags[p] = sf2;
  }

}

void PandaLeptonicAnalyzer::RegisterTrigger(TString path, std::vector<unsigned> &idxs) {
  unsigned idx = event.registerTrigger(path);
  if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::RegisterTrigger",
            TString::Format("At %u found trigger=%s",idx,path.Data()));
  idxs.push_back(idx);
}

// run
void PandaLeptonicAnalyzer::Run() {

  // INITIALIZE --------------------------------------------------------------------------
  unsigned int nEvents = tIn->GetEntries();
  unsigned int nZero = 0;
  if (lastEvent>=0 && lastEvent<(int)nEvents)
    nEvents = lastEvent;
  if (firstEvent>=0)
    nZero = firstEvent;

  if (!fOut || !tIn) {
    PError("PandaLeptonicAnalyzer::Run","NOT SETUP CORRECTLY");
    exit(1);
  }

  panda::JetCollection* jets(0);
  jets = &event.chsAK4Jets;

  // these are bins of b-tagging eff in pT and eta, derived in 8024 TT MC
  std::vector<double> vbtagpt {20.0,50.0,80.0,120.0,200.0,300.0,400.0,500.0,700.0,1000.0};
  std::vector<double> vbtageta {0.0,0.5,1.5,2.5};
  std::vector<std::vector<double>> lfeff  = {{0.081,0.065,0.060,0.063,0.072,0.085,0.104,0.127,0.162},
                       {0.116,0.097,0.092,0.099,0.112,0.138,0.166,0.185,0.222},
                       {0.173,0.145,0.149,0.175,0.195,0.225,0.229,0.233,0.250}};
  std::vector<std::vector<double>> ceff = {{0.377,0.389,0.391,0.390,0.391,0.375,0.372,0.392,0.435},
                      {0.398,0.407,0.416,0.424,0.424,0.428,0.448,0.466,0.500},
                      {0.375,0.389,0.400,0.425,0.437,0.459,0.481,0.534,0.488}};
  std::vector<std::vector<double>> beff = {{0.791,0.815,0.825,0.835,0.821,0.799,0.784,0.767,0.760},
                      {0.794,0.816,0.829,0.836,0.823,0.804,0.798,0.792,0.789},
                      {0.739,0.767,0.780,0.789,0.776,0.771,0.779,0.787,0.806}};
  Binner btagpt(vbtagpt);
  Binner btageta(vbtageta);

  JetCorrectionUncertainty *uncReader=0;
  JetCorrectionUncertainty *uncReaderAK4=0;
  FactorizedJetCorrector *scaleReaderAK4=0;

  std::vector<unsigned int> metTriggers;
  std::vector<unsigned int> phoTriggers;
  std::vector<unsigned int> muegTriggers;
  std::vector<unsigned int> mumuTriggers;
  std::vector<unsigned int> muTriggers;
  std::vector<unsigned int> egegTriggers;
  std::vector<unsigned int> egTriggers;

  if (1) {
    std::vector<TString> metTriggerPaths = {
          "HLT_PFMET170_NoiseCleaned",
          "HLT_PFMET170_HBHECleaned",
          "HLT_PFMET170_JetIdCleaned",
          "HLT_PFMET170_NotCleaned",
          "HLT_PFMET170_HBHE_BeamHaloCleaned",
          "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight",
          "HLT_PFMETNoMu110_NoiseCleaned_PFMHTNoMu110_IDTight",
          "HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight",
          "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight",
          "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight",
          "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
          "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"
    };
    std::vector<TString> phoTriggerPaths = {
          "HLT_Photon175",
          "HLT_Photon165_HE10",
          "HLT_Photon36_R9Id90_HE10_IsoM",
          "HLT_Photon50_R9Id90_HE10_IsoM",
          "HLT_Photon75_R9Id90_HE10_IsoM",
          "HLT_Photon90_R9Id90_HE10_IsoM",
          "HLT_Photon120_R9Id90_HE10_IsoM",
          "HLT_Photon165_R9Id90_HE10_IsoM",
          "HLT_Photon300_NoHE",

          "HLT_ECALHT800"
    };

    std::vector<TString> muegTriggerPaths = {
	  "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
	  "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
	  "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
	  "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
	  "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
	  "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
	  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
	  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"
    };

    std::vector<TString> mumuTriggerPaths = {
	  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
	  "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
	  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
	  "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"

    };

    std::vector<TString> muTriggerPaths = {
	  "HLT_IsoMu24",
	  "HLT_IsoTkMu24",
	  "HLT_IsoMu22",
	  "HLT_IsoTkMu22",
	  "HLT_Mu45_eta2p1",
	  "HLT_Mu50"

    };

    std::vector<TString> egegTriggerPaths = {
	  "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
	  "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf"
    };

    std::vector<TString> egTriggerPaths = {
          "HLT_Ele25_eta2p1_WPTight_Gsf",
	  "HLT_Ele27_eta2p1_WPLoose_Gsf",
	  "HLT_Ele27_WPTight_Gsf",
	  "HLT_Ele30_WPTight_Gsf",
	  "HLT_Ele35_WPLoose_Gsf",
          "HLT_Ele27_WP85_Gsf",
          "HLT_Ele27_WPLoose_Gsf",
          "HLT_Ele105_CaloIdVT_GsfTrkIdT",
          "HLT_Ele27_eta2p1_WPTight_Gsf",
          "HLT_Ele32_eta2p1_WPTight_Gsf",
          "HLT_ECALHT800"
    };

    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading MET triggers");
    for (auto path : metTriggerPaths) {
      RegisterTrigger(path,metTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading SinglePhoton triggers");
    for (auto path : phoTriggerPaths) {
      RegisterTrigger(path,phoTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading MuEG triggers");
    for (auto path : muegTriggerPaths) {
      RegisterTrigger(path,muegTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading MuMu triggers");
    for (auto path : mumuTriggerPaths) {
      RegisterTrigger(path,mumuTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading Mu triggers");
    for (auto path : muTriggerPaths) {
      RegisterTrigger(path,muTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading EGEG triggers");
    for (auto path : egegTriggerPaths) {
      RegisterTrigger(path,egegTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading EG triggers");
    for (auto path : egTriggerPaths) {
      RegisterTrigger(path,egTriggers);
    }

  }

  float EGMSCALE = isData ? 1 : 1;

  // set up reporters
  unsigned int iE=0;
  ProgressReporter pr("PandaLeptonicAnalyzer::Run",&iE,&nEvents,10);
  TimeReporter tr("PandaLeptonicAnalyzer::Run",DEBUG);

  bool applyJSON = flags["applyJSON"];

  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    tr.Start();
    pr.Report();
    ResetBranches();
    event.getEntry(*tIn,iE);

    tr.TriggerEvent(TString::Format("GetEntry %u",iE));
    if (DEBUG>2) {
      PDebug("PandaLeptonicAnalyzer::Run::Dump","");
      event.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaLeptonicAnalyzer::Run::Dump","");
      event.photons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaLeptonicAnalyzer::Run::Dump","");
      event.muons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaLeptonicAnalyzer::Run::Dump","");
      event.electrons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaLeptonicAnalyzer::Run::Dump","");
      event.chsAK4Jets.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaLeptonicAnalyzer::Run::Dump","");
      event.pfMet.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaLeptonicAnalyzer::Run::Dump","");
      event.metMuOnlyFix.print(std::cout, 2);
      std::cout << std::endl;
    }

    // event info
    gt->runNumber = event.runNumber;
    gt->lumiNumber = event.lumiNumber;
    gt->eventNumber = event.eventNumber;
    gt->npv = event.npv;
    gt->pu = event.npvTrue;
    gt->mcWeight = event.weight;
    gt->metFilter = (event.metFilters.pass()) ? 1 : 0;
    // these two are not need since we use muon-fixed MET
    // gt->metFilter = (gt->metFilter==1 && !event.metFilters.badMuons) ? 1 : 0;
    // gt->metFilter = (gt->metFilter==1 && !event.metFilters.duplicateMuons) ? 1 : 0;
    gt->metFilter = (gt->metFilter==1 && !event.metFilters.badPFMuons) ? 1 : 0;
    gt->metFilter = (gt->metFilter==1 && !event.metFilters.badChargedHadrons) ? 1 : 0;
    gt->egmFilter = (!event.metFilters.dupECALClusters) ? 1 : 0;
    gt->egmFilter = (gt->egmFilter==1 && !event.metFilters.unfixedECALHits) ? 1 : 0;

    // save triggers
    for (auto iT : metTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kMETTrig;
      break;
     }
    }
    for (auto iT : phoTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kSinglePhoTrig;
      break;
     }
    }
    for (auto iT : muegTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kMuEGTrig;
      break;
     }
    }
    for (auto iT : mumuTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kMuMuTrig;
      break;
     }
    }
    for (auto iT : muTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kMuTrig;
      break;
     }
    }
    for (auto iT : egegTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kEGEGTrig;
      break;
     }
    }
    for (auto iT : egTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kEGTrig;
      break;
     }
    }

    if (isData) {
      // check the json
      if (applyJSON && !PassGoodLumis(gt->runNumber,gt->lumiNumber))
        continue;

    } else {
      gt->sf_pu     = GetCorr(cPU    ,gt->pu);
      gt->sf_puUp   = GetCorr(cPUUp  ,gt->pu);
      gt->sf_puDown = GetCorr(cPUDown,gt->pu);
    }

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

    tr.TriggerEvent("initialize");

    // met
    gt->pfmet = event.pfMet.pt;
    gt->pfmetphi = event.pfMet.phi;
    gt->pfmetRaw = event.rawMet.pt;
    gt->pfmetUp = event.pfMet.ptCorrUp;
    gt->pfmetDown = event.pfMet.ptCorrDown;
    gt->puppimet = event.puppiMet.pt;
    gt->puppimetphi = event.puppiMet.phi;
    gt->calomet = event.caloMet.pt;
    gt->calometphi = event.caloMet.phi;
    gt->trkmet = event.trkMet.pt;
    gt->trkmetphi = event.trkMet.phi;
    TLorentzVector vPFMET, vPuppiMET;
    vPFMET.SetPtEtaPhiM(gt->pfmet,0,gt->pfmetphi,0);
    vPuppiMET.SetPtEtaPhiM(gt->puppimet,0,gt->puppimetphi,0);
    TVector2 vMETNoMu; vMETNoMu.SetMagPhi(gt->pfmet,gt->pfmetphi); //       for trigger eff

    tr.TriggerEvent("met");

    //electrons
    std::vector<panda::Lepton*> looseLeps;
    for (auto& ele : event.electrons) {
     float pt = ele.pt()*EGMSCALE; float eta = ele.eta(); float aeta = fabs(eta);
      if (pt<10 || aeta>2.5)
        continue;
      if (!ele.veto)
        continue;
      if (!ElectronIP(ele.eta(),ele.dxy,ele.dz))
        continue;
      looseLeps.push_back(&ele);
      matchLeps.push_back(&ele);
    }

    // muons
    for (auto& mu : event.muons) {
     float pt = mu.pt(); float eta = mu.eta(); float aeta = fabs(eta);
      if (pt<10 || aeta>2.4)
        continue;
      if (!mu.loose)
        continue;
      looseLeps.push_back(&mu);
      matchLeps.push_back(&mu);
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
     int nToSort = gt->nLooseLep;
     std::partial_sort(looseLeps.begin(),looseLeps.begin()+nToSort,looseLeps.end(),ptsort);
    }
    int lep_counter=1;
    for (auto* lep : looseLeps) {
      if      (lep_counter==1) {
       gt->looseLep1Pt  = lep->pt();
       gt->looseLep1Eta = lep->eta();
       gt->looseLep1Phi = lep->phi();
      }
      else if (lep_counter==2) {
       gt->looseLep2Pt  = lep->pt();
       gt->looseLep2Eta = lep->eta();
       gt->looseLep2Phi = lep->phi();
      } 
      else if (lep_counter==3) {
       gt->looseLep3Pt  = lep->pt();
       gt->looseLep3Eta = lep->eta();
       gt->looseLep3Phi = lep->phi();
      } 
      else if (lep_counter==4) {
       gt->looseLep4Pt  = lep->pt();
       gt->looseLep4Eta = lep->eta();
       gt->looseLep4Phi = lep->phi();
      } 
      else {
        break;
      }
      // now specialize lepton types
      panda::Muon *mu = dynamic_cast<panda::Muon*>(lep);
      if (mu!=NULL) {
        bool isFake   = mu->tight  && mu->combIso()/mu->pt() < 0.4 && mu->chIso/mu->pt() < 0.4;
        bool isMedium = mu->medium && mu->combIso()/mu->pt() < 0.15;
        bool isTight  = mu->tight  && mu->combIso()/mu->pt() < 0.15;
        if      (lep_counter==1) {
          gt->looseLep1PdgId = mu->charge*-13;
                       gt->looseLep1SelBit |= kLoose;
          if(isFake)   gt->looseLep1SelBit |= kFake;
          if(isMedium) gt->looseLep1SelBit |= kMedium;
          if(isTight)  gt->looseLep1SelBit |= kTight;
	  gt->sf_trk1    = GetCorr(cTrackingMuon,mu->eta());
	  gt->sf_loose1  = GetCorr(cLooseMuonId,TMath::Abs(mu->eta()),mu->pt())  * GetCorr(cLooseMuonIso,TMath::Abs(mu->eta()),mu->pt());
	  gt->sf_medium1 = GetCorr(cMediumMuonId,TMath::Abs(mu->eta()),mu->pt()) * GetCorr(cMediumMuonIso,TMath::Abs(mu->eta()),mu->pt());
	  gt->sf_tight1  = GetCorr(cTightMuonId,TMath::Abs(mu->eta()),mu->pt())  * GetCorr(cTightMuonIso,TMath::Abs(mu->eta()),mu->pt());
        }
	else if (lep_counter==2) {
          gt->looseLep2PdgId = mu->charge*-13;
                       gt->looseLep2SelBit |= kLoose;
          if(isFake)   gt->looseLep2SelBit |= kFake;
          if(isMedium) gt->looseLep2SelBit |= kMedium;
          if(isTight)  gt->looseLep2SelBit |= kTight;
	  gt->sf_trk2    = GetCorr(cTrackingMuon,mu->eta());
	  gt->sf_loose2  = GetCorr(cLooseMuonId,TMath::Abs(mu->eta()),mu->pt())  * GetCorr(cLooseMuonIso,TMath::Abs(mu->eta()),mu->pt());
	  gt->sf_medium2 = GetCorr(cMediumMuonId,TMath::Abs(mu->eta()),mu->pt()) * GetCorr(cMediumMuonIso,TMath::Abs(mu->eta()),mu->pt());
	  gt->sf_tight2  = GetCorr(cTightMuonId,TMath::Abs(mu->eta()),mu->pt())  * GetCorr(cTightMuonIso,TMath::Abs(mu->eta()),mu->pt());
        }
	else if (lep_counter==3) {
          gt->looseLep3PdgId = mu->charge*-13;
                       gt->looseLep3SelBit |= kLoose;
          if(isFake)   gt->looseLep3SelBit |= kFake;
          if(isMedium) gt->looseLep3SelBit |= kMedium;
          if(isTight)  gt->looseLep3SelBit |= kTight;
	  gt->sf_trk3    = GetCorr(cTrackingMuon,mu->eta());
	  gt->sf_loose3  = GetCorr(cLooseMuonId,TMath::Abs(mu->eta()),mu->pt())  * GetCorr(cLooseMuonIso,TMath::Abs(mu->eta()),mu->pt());
	  gt->sf_medium3 = GetCorr(cMediumMuonId,TMath::Abs(mu->eta()),mu->pt()) * GetCorr(cMediumMuonIso,TMath::Abs(mu->eta()),mu->pt());
	  gt->sf_tight3  = GetCorr(cTightMuonId,TMath::Abs(mu->eta()),mu->pt())  * GetCorr(cTightMuonIso,TMath::Abs(mu->eta()),mu->pt());
        }
	else if (lep_counter==4) {
          gt->looseLep4PdgId = mu->charge*-13;
                       gt->looseLep4SelBit |= kLoose;
          if(isFake)   gt->looseLep4SelBit |= kFake;
          if(isMedium) gt->looseLep4SelBit |= kMedium;
          if(isTight)  gt->looseLep4SelBit |= kTight;
	  gt->sf_trk4    = GetCorr(cTrackingMuon,mu->eta());
	  gt->sf_loose4  = GetCorr(cLooseMuonId,TMath::Abs(mu->eta()),mu->pt())  * GetCorr(cLooseMuonIso,TMath::Abs(mu->eta()),mu->pt());
	  gt->sf_medium4 = GetCorr(cMediumMuonId,TMath::Abs(mu->eta()),mu->pt()) * GetCorr(cMediumMuonIso,TMath::Abs(mu->eta()),mu->pt());
	  gt->sf_tight4  = GetCorr(cTightMuonId,TMath::Abs(mu->eta()),mu->pt())  * GetCorr(cTightMuonIso,TMath::Abs(mu->eta()),mu->pt());
        }
      } else {
        panda::Electron *ele = dynamic_cast<panda::Electron*>(lep);
        bool isFake   = ele->hltsafe;
        bool isMedium = ele->medium;
        bool isTight  = ele->tight;
        if      (lep_counter==1) {
          gt->looseLep1Pt *= EGMSCALE;
          gt->looseLep1PdgId = ele->charge*-11;
                       gt->looseLep1SelBit |= kLoose;
          if(isFake)   gt->looseLep1SelBit |= kFake;
          if(isMedium) gt->looseLep1SelBit |= kMedium;
          if(isTight)  gt->looseLep1SelBit |= kTight;
	  gt->sf_trk1    = GetCorr(cTrackingElectron,ele->eta(),ele->pt());
	  gt->sf_loose1  = GetCorr(cLooseElectronId,ele->eta(),ele->pt());
	  gt->sf_medium1 = GetCorr(cMediumElectronId,ele->eta(),ele->pt());
	  gt->sf_tight1  = GetCorr(cTightElectronId,ele->eta(),ele->pt());
        } 
	else if (lep_counter==2) {
          gt->looseLep2Pt *= EGMSCALE;
          gt->looseLep2PdgId = ele->charge*-11;
                       gt->looseLep2SelBit |= kLoose;
          if(isFake)   gt->looseLep2SelBit |= kFake;
          if(isMedium) gt->looseLep2SelBit |= kMedium;
          if(isTight)  gt->looseLep2SelBit |= kTight;
	  gt->sf_trk2    = GetCorr(cTrackingElectron,ele->eta(),ele->pt());
	  gt->sf_loose2  = GetCorr(cLooseElectronId,ele->eta(),ele->pt());
	  gt->sf_medium2 = GetCorr(cMediumElectronId,ele->eta(),ele->pt());
	  gt->sf_tight2  = GetCorr(cTightElectronId,ele->eta(),ele->pt());
        }
	else if (lep_counter==3) {
          gt->looseLep3Pt *= EGMSCALE;
          gt->looseLep3PdgId = ele->charge*-11;
                       gt->looseLep3SelBit |= kLoose;
          if(isFake)   gt->looseLep3SelBit |= kFake;
          if(isMedium) gt->looseLep3SelBit |= kMedium;
          if(isTight)  gt->looseLep3SelBit |= kTight;
	  gt->sf_trk3    = GetCorr(cTrackingElectron,ele->eta(),ele->pt());
	  gt->sf_loose3  = GetCorr(cLooseElectronId,ele->eta(),ele->pt());
	  gt->sf_medium3 = GetCorr(cMediumElectronId,ele->eta(),ele->pt());
	  gt->sf_tight3  = GetCorr(cTightElectronId,ele->eta(),ele->pt());
        }
	else if (lep_counter==4) {
          gt->looseLep4Pt *= EGMSCALE;
          gt->looseLep4PdgId = ele->charge*-11;
                       gt->looseLep4SelBit |= kLoose;
          if(isFake)   gt->looseLep4SelBit |= kFake;
          if(isMedium) gt->looseLep4SelBit |= kMedium;
          if(isTight)  gt->looseLep4SelBit |= kTight;
	  gt->sf_trk4    = GetCorr(cTrackingElectron,ele->eta(),ele->pt());
	  gt->sf_loose4  = GetCorr(cLooseElectronId,ele->eta(),ele->pt());
	  gt->sf_medium4 = GetCorr(cMediumElectronId,ele->eta(),ele->pt());
	  gt->sf_tight4  = GetCorr(cTightElectronId,ele->eta(),ele->pt());
        }
      }
      ++lep_counter;
    }

    tr.TriggerEvent("leptons");

    // photons
    gt->nLoosePhoton = 0;
    for (auto& pho : event.photons) {
      if (!pho.medium || !pho.pixelVeto)
        continue;
      float pt = pho.pt() * EGMSCALE;
      if (pt<1) continue;
      float eta = pho.eta(), phi = pho.phi();
      if (pt<20 || fabs(eta)>2.5)
        continue;
      if (IsMatched(&matchLeps,0.16,eta,phi))
        continue;
      gt->nLoosePhoton++;
      matchPhos.push_back(&pho);
      if (gt->nLoosePhoton==1) {
        gt->loosePho1Pt = pt;
        gt->loosePho1Eta = eta;
        gt->loosePho1Phi = phi;
      }
    }

    tr.TriggerEvent("photons");

    // first identify interesting jets
    vector<panda::Jet*> cleanedJets;
    vector<int> btagindices;
    TLorentzVector vJet;
    panda::Jet *jet1=0, *jet2=0, *jet3=0, *jet4=0;
    panda::Jet *jetUp1=0, *jetUp2=0, *jetUp3=0, *jetUp4=0;
    panda::Jet *jetDown1=0, *jetDown2=0, *jetDown3=0, *jetDown4=0;
    gt->dphipuppimet=999; gt->dphipfmet=999;
    float maxJetEta = 4.7;
    unsigned nJetDPhi = 1;

    for (auto& jet : *jets) {

      // only do eta-phi checks here
      if (abs(jet.eta()) > maxJetEta)
         continue;
      if (IsMatched(&matchLeps,0.16,jet.eta(),jet.phi()))
         continue;

      bool isLoose = jet.loose;
      bool isTight = jet.tight;

      if (jet.pt()>20 && jet.csv>0.5426) ++(gt->jetNLBtags);
      if (jet.pt()>20 && jet.csv>0.8484) ++(gt->jetNMBtags);
      if (jet.pt()>20 && jet.csv>0.9535) ++(gt->jetNLBtags);

      if (jet.pt()>30) { // nominal jets
	cleanedJets.push_back(&jet);
	if      (cleanedJets.size()==1) {
          jet1 = &jet;
          gt->jet1Pt   = jet.pt();
          gt->jet1Eta  = jet.eta();
          gt->jet1Phi  = jet.phi();
          gt->jet1BTag = jet.csv;
          if(isLoose) gt->jet1SelBit |= kLoose;
          if(isTight) gt->jet1SelBit |= kTight;
	}
	else if (cleanedJets.size()==2) {
          jet2 = &jet;
          gt->jet2Pt   = jet.pt();
          gt->jet2Eta  = jet.eta();
          gt->jet2Phi  = jet.phi();
          gt->jet2BTag = jet.csv;
          if(isLoose) gt->jet2SelBit |= kLoose;
          if(isTight) gt->jet2SelBit |= kTight;
	}
	else if (cleanedJets.size()==3) {
          jet3 = &jet;
          gt->jet3Pt   = jet.pt();
          gt->jet3Eta  = jet.eta();
          gt->jet3Phi  = jet.phi();
          gt->jet3BTag = jet.csv;
          if(isLoose) gt->jet3SelBit |= kLoose;
          if(isTight) gt->jet3SelBit |= kTight;
	}
	else if (cleanedJets.size()==4) {
          jet4 = &jet;
          gt->jet4Pt   = jet.pt();
          gt->jet4Eta  = jet.eta();
          gt->jet4Phi  = jet.phi();
          gt->jet4BTag = jet.csv;
          if(isLoose) gt->jet4SelBit |= kLoose;
          if(isTight) gt->jet4SelBit |= kTight;
	}

	// compute dphi wrt mets
	if (cleanedJets.size() <= nJetDPhi) {
          vJet.SetPtEtaPhiM(jet.pt(),jet.eta(),jet.phi(),jet.m());
          gt->dphipuppimet = std::min(fabs(vJet.DeltaPhi(vPuppiMET)),(double)gt->dphipuppimet);
          gt->dphipfmet = std::min(fabs(vJet.DeltaPhi(vPFMET)),(double)gt->dphipfmet);
	}
      }

      // do jes variation OUTSIDE of pt>30 check
      if (jet.ptCorrUp>30) {
	if (jet.ptCorrUp > gt->jet1PtUp) {
          if (jetUp1) {
            jetUp4 = jetUp3;
            gt->jet4PtUp  = gt->jet3PtUp;
            gt->jet4EtaUp = gt->jet3EtaUp;
            jetUp3 = jetUp2;
            gt->jet3PtUp  = gt->jet2PtUp;
            gt->jet3EtaUp = gt->jet2EtaUp;
            jetUp2 = jetUp1;
            gt->jet2PtUp  = gt->jet1PtUp;
            gt->jet2EtaUp = gt->jet1EtaUp;
          }
          jetUp1 = &jet;
          gt->jet1PtUp  = jet.ptCorrUp;
          gt->jet1EtaUp = jet.eta();
	}
	else if (jet.ptCorrUp > gt->jet2PtUp) {
          if (jetUp2) {
            jetUp4 = jetUp3;
            gt->jet4PtUp  = gt->jet3PtUp;
            gt->jet4EtaUp = gt->jet3EtaUp;
            jetUp3 = jetUp2;
            gt->jet3PtUp  = gt->jet2PtUp;
            gt->jet3EtaUp = gt->jet2EtaUp;
          }
          jetUp2 = &jet;
          gt->jet2PtUp  = jet.ptCorrUp;
          gt->jet2EtaUp = jet.eta();
	}
	else if (jet.ptCorrUp > gt->jet3PtUp) {
          if (jetUp3) {
            jetUp4 = jetUp3;
            gt->jet4PtUp  = gt->jet3PtUp;
            gt->jet4EtaUp = gt->jet3EtaUp;
          }
          jetUp3 = &jet;
          gt->jet3PtUp  = jet.ptCorrUp;
          gt->jet3EtaUp = jet.eta();
	}
	else if (jet.ptCorrUp > gt->jet4PtUp) {
          jetUp4 = &jet;
          gt->jet4PtUp  = jet.ptCorrUp;
          gt->jet4EtaUp = jet.eta();
	}
      }

      if (jet.ptCorrDown>30) {
	if (jet.ptCorrDown > gt->jet1PtDown) {
          if (jetDown1) {
            jetDown4 = jetDown3;
            gt->jet4PtDown  = gt->jet3PtDown;
            gt->jet4EtaDown = gt->jet3EtaDown;
            jetDown3 = jetDown2;
            gt->jet3PtDown  = gt->jet2PtDown;
            gt->jet3EtaDown = gt->jet2EtaDown;
            jetDown2 = jetDown1;
            gt->jet2PtDown  = gt->jet1PtDown;
            gt->jet2EtaDown = gt->jet1EtaDown;
          }
          jetDown1 = &jet;
          gt->jet1PtDown  = jet.ptCorrDown;
          gt->jet1EtaDown = jet.eta();
	}
	else if (jet.ptCorrDown > gt->jet2PtDown) {
          if (jetDown2) {
            jetDown4 = jetDown3;
            gt->jet4PtDown  = gt->jet3PtDown;
            gt->jet4EtaDown = gt->jet3EtaDown;
            jetDown3 = jetDown2;
            gt->jet3PtDown  = gt->jet2PtDown;
            gt->jet3EtaDown = gt->jet2EtaDown;
          }
          jetDown2 = &jet;
          gt->jet2PtDown  = jet.ptCorrDown;
          gt->jet2EtaDown = jet.eta();
	}
	else if (jet.ptCorrDown > gt->jet3PtDown) {
          if (jetDown3) {
            jetDown4 = jetDown3;
            gt->jet4PtDown  = gt->jet3PtDown;
            gt->jet4EtaDown = gt->jet3EtaDown;
          }
          jetDown3 = &jet;
          gt->jet3PtDown  = jet.ptCorrDown;
          gt->jet3EtaDown = jet.eta();
	}
	else if (jet.ptCorrDown > gt->jet4PtDown) {
          jetDown4 = &jet;
          gt->jet4PtDown  = jet.ptCorrDown;
          gt->jet4EtaDown = jet.eta();
	}
      }
    } // Jet loop

    gt->nJet = cleanedJets.size();

    tr.TriggerEvent("jets");

    for (auto& tau : event.taus) {
      if (!tau.decayMode || !tau.decayModeNew)
        continue;
      if (!tau.looseIsoMVA)
        continue;
      if (tau.pt()<18 || fabs(tau.eta())>2.3)
        continue;
      if (IsMatched(&matchLeps,0.16,tau.eta(),tau.phi()))
        continue;
      gt->nTau++;
    }

    tr.TriggerEvent("taus");

    gt->sf_tt = 1;
    gt->genLep1Pt = 0;
    gt->genLep1Eta = -1;
    gt->genLep1Phi = -1;
    gt->genLep1PdgId = 0;
    gt->genLep2Pt = 0;
    gt->genLep2Eta = -1;
    gt->genLep2Phi = -1;
    gt->genLep2PdgId = 0;
    gt->looseGenLep1PdgId = 0;
    gt->looseGenLep2PdgId = 0;
    gt->looseGenLep3PdgId = 0;
    gt->looseGenLep4PdgId = 0;
    TLorentzVector v1,v2,v3,v4;
    if (gt->nLooseLep>=1) {
      panda::Lepton *lep1=looseLeps[0];
      v1.SetPtEtaPhiM(lep1->pt(),lep1->eta(),lep1->phi(),lep1->m());
    }
    if (gt->nLooseLep>=2) {
      panda::Lepton *lep2=looseLeps[1];
      v2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());
    }
    if (gt->nLooseLep>=3) {
      panda::Lepton *lep3=looseLeps[2];
      v3.SetPtEtaPhiM(lep3->pt(),lep3->eta(),lep3->phi(),lep3->m());
    }
    if (gt->nLooseLep>=4) {
      panda::Lepton *lep4=looseLeps[3];
      v4.SetPtEtaPhiM(lep4->pt(),lep4->eta(),lep4->phi(),lep4->m());
    }
    // gen lepton matching
    if (!isData) {
      std::vector<int> targetsLepton;
      std::vector<int> targetsPhoton;
      std::vector<int> targetsV;
      std::vector<int> targetsTop;
      std::vector<int> targetsN;

      int nGen = event.genParticles.size();
      for (int iG=0; iG!=nGen; ++iG) {
        auto& part(event.genParticles.at(iG));
        int pdgid = part.pdgid;
        unsigned int abspdgid = abs(pdgid);
        if ((abspdgid == 11 || abspdgid == 13) &&
	    (part.statusFlags == GenParticle::kIsPrompt || 
	     part.statusFlags == GenParticle::kIsTauDecayProduct || part.statusFlags == GenParticle::kIsPromptTauDecayProduct || 
	     part.statusFlags == GenParticle::kIsDirectTauDecayProduct || part.statusFlags == GenParticle::kIsDirectPromptTauDecayProduct ))
          targetsLepton.push_back(iG);

        if (abspdgid == 22)
          targetsPhoton.push_back(iG);

        if (abspdgid == 23 || abspdgid == 24)
          targetsV.push_back(iG);

        if (abspdgid == 6)
          targetsTop.push_back(iG);

        if (abspdgid == 12 || abspdgid == 14 || abspdgid == 16)
          targetsN.push_back(iG);

      } //looking for targets

      TLorentzVector the_rhoP4(0,0,0,0);
      double bosonPtMin = 1000000000;
      for (int iG : targetsLepton) {
        auto& part(event.genParticles.at(iG));
        TLorentzVector dressedLepton;
        dressedLepton.SetPtEtaPhiM(part.pt(),part.eta(),part.phi(),part.m());

        // check there is no further copy:
        bool isLastCopy=true;
        for (int kG : targetsLepton) {
          if (event.genParticles.at(kG).parent.get() == &part) {
            isLastCopy=false;
            break;
          }
        }
        if (!isLastCopy)
          continue;
	
	the_rhoP4 = the_rhoP4 + dressedLepton;
	
        for (int jG : targetsPhoton) {
          auto& partj(event.genParticles.at(jG));

          // check there is no further copy:
          bool isLastCopy=true;
          for (int kG : targetsPhoton) {
            if (event.genParticles.at(kG).parent.get() == &part) {
              isLastCopy=false;
              break;
            }
          }
          if (!isLastCopy)
            continue;
	  
	  if(abs(partj.pdgid) == 22 && DeltaR2(part.eta(),part.phi(),partj.eta(),partj.phi()) < 0.1*0.1) {
            TLorentzVector photonV;
            photonV.SetPtEtaPhiM(partj.pt(),partj.eta(),partj.phi(),partj.m());
	    dressedLepton += photonV;
	  }
	}

	if     (dressedLepton.Pt() > gt->genLep1Pt){
          gt->genLep2Pt    = gt->genLep1Pt; 
          gt->genLep2Eta   = gt->genLep1Eta;
          gt->genLep2Phi   = gt->genLep1Phi;
          gt->genLep2PdgId = gt->genLep1PdgId; 
          gt->genLep1Pt    = dressedLepton.Pt();
          gt->genLep1Eta   = dressedLepton.Eta();
          gt->genLep1Phi   = dressedLepton.Phi();
          gt->genLep1PdgId = part.pdgid;
        }
	else if(dressedLepton.Pt() > gt->genLep2Pt){
          gt->genLep2Pt    = dressedLepton.Pt();
          gt->genLep2Eta   = dressedLepton.Eta();
          gt->genLep2Phi   = dressedLepton.Phi();
          gt->genLep2PdgId = part.pdgid; 
        }
 
        if(v1.Pt() > 0 && DeltaR2(part.eta(),part.phi(),v1.Eta(),v1.Phi()) < 0.1*0.1) {
	  if     (part.statusFlags == GenParticle::kIsTauDecayProduct || part.statusFlags == GenParticle::kIsPromptTauDecayProduct || 
	          part.statusFlags == GenParticle::kIsDirectTauDecayProduct || part.statusFlags == GenParticle::kIsDirectPromptTauDecayProduct) gt->looseGenLep1PdgId = 2;
	  else if(part.statusFlags == GenParticle::kIsPrompt) gt->looseGenLep1PdgId = 1;
	  if(part.pdgid != gt->looseLep1PdgId) gt->looseGenLep1PdgId = -1 * gt->looseGenLep1PdgId;
	}

        if(v2.Pt() > 0 && DeltaR2(part.eta(),part.phi(),v2.Eta(),v2.Phi()) < 0.1*0.1) {
	  if     (part.statusFlags == GenParticle::kIsTauDecayProduct || part.statusFlags == GenParticle::kIsPromptTauDecayProduct || 
	          part.statusFlags == GenParticle::kIsDirectTauDecayProduct || part.statusFlags == GenParticle::kIsDirectPromptTauDecayProduct) gt->looseGenLep2PdgId = 2;
	  else if(part.statusFlags == GenParticle::kIsPrompt) gt->looseGenLep2PdgId = 1;
	  if(part.pdgid != gt->looseLep2PdgId) gt->looseGenLep2PdgId = -1 * gt->looseGenLep2PdgId;
	}

        if(v3.Pt() > 0 && DeltaR2(part.eta(),part.phi(),v3.Eta(),v3.Phi()) < 0.1*0.1) {
	  if     (part.statusFlags == GenParticle::kIsTauDecayProduct || part.statusFlags == GenParticle::kIsPromptTauDecayProduct || 
	          part.statusFlags == GenParticle::kIsDirectTauDecayProduct || part.statusFlags == GenParticle::kIsDirectPromptTauDecayProduct) gt->looseGenLep3PdgId = 2;
	  else if(part.statusFlags == GenParticle::kIsPrompt) gt->looseGenLep3PdgId = 1;
	  if(part.pdgid != gt->looseLep3PdgId) gt->looseGenLep3PdgId = -1 * gt->looseGenLep3PdgId;
	}

        if(v4.Pt() > 0 && DeltaR2(part.eta(),part.phi(),v4.Eta(),v4.Phi()) < 0.1*0.1) {
	  if     (part.statusFlags == GenParticle::kIsTauDecayProduct || part.statusFlags == GenParticle::kIsPromptTauDecayProduct || 
	          part.statusFlags == GenParticle::kIsDirectTauDecayProduct || part.statusFlags == GenParticle::kIsDirectPromptTauDecayProduct) gt->looseGenLep4PdgId = 2;
	  else if(part.statusFlags == GenParticle::kIsPrompt) gt->looseGenLep4PdgId = 1;
	  if(part.pdgid != gt->looseLep4PdgId) gt->looseGenLep4PdgId = -1 * gt->looseGenLep4PdgId;
	}
      }
      
      // Filling dilepton Pt at gen level
      if(gt->genLep1Pt > 25 && TMath::Abs(gt->genLep1Eta) < 2.5 && 
         gt->genLep2Pt > 25 && TMath::Abs(gt->genLep2Eta) < 2.5){
        TLorentzVector genlep1;
        genlep1.SetPtEtaPhiM(gt->genLep1Pt,gt->genLep1Eta,gt->genLep1Phi,0.0);
        TLorentzVector genlep2;
        genlep2.SetPtEtaPhiM(gt->genLep2Pt,gt->genLep2Eta,gt->genLep2Phi,0.0);
	TLorentzVector dilep = genlep1 + genlep2;
	if(TMath::Abs(dilep.M()-91.1876) < 15.0) {
	  if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtMM->Fill(dilep.Pt(),event.weight);
	  else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtEE->Fill(dilep.Pt(),event.weight);
	  if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13 && dilep.Pt() < 50.0) hDDilLowPtMM->Fill(dilep.Pt(),event.weight);
	  else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11 && dilep.Pt() < 50.0) hDDilLowPtEE->Fill(dilep.Pt(),event.weight);
	  if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPt2MM->Fill(dilep.Pt()*dilep.Pt(),event.weight);
	  else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPt2EE->Fill(dilep.Pt()*dilep.Pt(),event.weight);
	}
      }

      for (int iG : targetsN) {
        auto& part(event.genParticles.at(iG));
        TLorentzVector neutrino;
        neutrino.SetPtEtaPhiM(part.pt(),part.eta(),part.phi(),part.m());

        // check there is no further copy:
        bool isLastCopy=true;
        for (int kG : targetsN) {
          if (event.genParticles.at(kG).parent.get() == &part) {
            isLastCopy=false;
            break;
          }
        }
        if (!isLastCopy)
          continue;
	
	the_rhoP4 = the_rhoP4 + neutrino;
      }

      TLorentzVector theZBosons(0,0,0,0);
      TLorentzVector theWBosons(0,0,0,0);
      int nZBosons = 0; int nWBosons = 0;
      for (int iG : targetsV) {
        auto& part(event.genParticles.at(iG));
        TLorentzVector boson;
        boson.SetPtEtaPhiM(part.pt(),part.eta(),part.phi(),part.m());

        // check there is no further copy:
        bool isLastCopy=true;
        for (int kG : targetsV) {
          if (event.genParticles.at(kG).parent.get() == &part) {
            isLastCopy=false;
            break;
          }
        }
        if (!isLastCopy)
          continue;
	
	if(boson.Pt() < bosonPtMin) bosonPtMin = boson.Pt();
	
	if(abs(part.pdgid) == 23) {theZBosons = theZBosons + boson; nZBosons++;}
	if(abs(part.pdgid) == 24) {theWBosons = theWBosons + boson; nWBosons++;}
      }
      if(nZBosons+nWBosons == 0) bosonPtMin = 0;

      if(nZBosons >= 2) {
        double the_rho = 0.0; if(the_rhoP4.P() > 0) the_rho = the_rhoP4.Pt()/the_rhoP4.P();
        double theZZCorr[2] {1,1};
        theZZCorr[0] = weightEWKCorr(bosonPtMin,1);
        float GENmZZ = theZBosons.M();
        theZZCorr[1] = kfactor_qqZZ_qcd_M(GENmZZ);
        gt->sf_zz = theZZCorr[0]*theZZCorr[1];
        if(the_rho <= 0.3) gt->sf_zzUnc = (1.0+TMath::Abs((theZZCorr[0]-1)*(15.99/9.89-1)));
	else               gt->sf_zzUnc = (1.0+TMath::Abs((theZZCorr[0]-1)               ));
      } else {
        gt->sf_zz    = 1.0;
	gt->sf_zzUnc = 1.0;
      }

      if(nWBosons == 1 && nZBosons == 1) {
        TLorentzVector theWZBoson = theWBosons + theZBosons;
        gt->sf_wz = weightEWKWZCorr(theWZBoson.M());
      } else {
        gt->sf_wz = 1.0;
      }
      
      if(nZBosons == 1) {
        gt->sf_zh     = weightZHEWKCorr(h1Corrs[ZHEwkCorr],     theZBosons.Pt());
        gt->sf_zhUp   = weightZHEWKCorr(h1Corrs[ZHEwkCorrUp],   theZBosons.Pt());
        gt->sf_zhDown = weightZHEWKCorr(h1Corrs[ZHEwkCorrDown], theZBosons.Pt());
      }
      else {
        gt->sf_zh     = 1.0;
        gt->sf_zhUp   = 1.0;
        gt->sf_zhDown = 1.0;
      }

      // ttbar pT weight
      TLorentzVector vT,vTbar;
      float pt_t=0, pt_tbar=0;
      for (int iG : targetsTop) {
        auto& part(event.genParticles.at(iG));

        // check there is no further copy:
        bool isLastCopy=true;
        for (int kG : targetsTop) {
          if (event.genParticles.at(kG).parent.get() == &part) {
            isLastCopy=false;
            break;
          }
        }
        if (!isLastCopy)
          continue;

        if (part.pdgid>0) {
         pt_t = part.pt();
        } else {
         pt_tbar = part.pt();
        }
      }
      if (pt_t>0 && pt_tbar>0) {
        gt->sf_tt = TMath::Sqrt(TMath::Exp(0.0615-0.0005*TMath::Min((float)400.,pt_t)) *
                  		TMath::Exp(0.0615-0.0005*TMath::Min((float)400.,pt_tbar)));
      }

    } // end gen study

    if (!PassPreselection())
      continue;

    tr.TriggerEvent("presel");

    if (!isData) {
      // now get the jet btag SFs
      vector<btagcand> btagcands;
      vector<btagcand> btagcands_alt;
      vector<double> sf_cent, sf_bUp, sf_bDown, sf_mUp, sf_mDown;
      vector<double> sf_cent_alt, sf_bUp_alt, sf_bDown_alt, sf_mUp_alt, sf_mDown_alt;

      unsigned int nJ = cleanedJets.size();
      for (unsigned int iJ=0; iJ!=nJ; ++iJ) {
        panda::Jet *jet = cleanedJets.at(iJ);
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
        if      (jet==cleanedJets.at(0)) {
          gt->jet1Flav = flavor;
          gt->jet1GenPt = genpt;
        } 
	else if (jet==cleanedJets.at(1)) {
          gt->jet2Flav = flavor;
          gt->jet2GenPt = genpt;
        }
	else if (jet==cleanedJets.at(2)) {
          gt->jet3Flav = flavor;
          gt->jet3GenPt = genpt;
        }
	else if (jet==cleanedJets.at(3)) {
          gt->jet4Flav = flavor;
          gt->jet4GenPt = genpt;
        }

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

      } // loop over jets

      EvalBTagSF(btagcands,sf_cent,GeneralLeptonicTree::bCent,GeneralLeptonicTree::bJet);
      EvalBTagSF(btagcands,sf_bUp,GeneralLeptonicTree::bBUp,GeneralLeptonicTree::bJet);
      EvalBTagSF(btagcands,sf_bDown,GeneralLeptonicTree::bBDown,GeneralLeptonicTree::bJet);
      EvalBTagSF(btagcands,sf_mUp,GeneralLeptonicTree::bMUp,GeneralLeptonicTree::bJet);
      EvalBTagSF(btagcands,sf_mDown,GeneralLeptonicTree::bMDown,GeneralLeptonicTree::bJet);

    }

    tr.TriggerEvent("ak4 gen-matching");

    // scale and PDF weights, if they exist
    gt->pdfUp = 1; gt->pdfDown = 1;
    if (!isData) {
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
      }
      tr.TriggerEvent("qcd uncertainties");

      unsigned nW = wIDs.size();
      if (nW) {
        for (unsigned iW=0; iW!=nW; ++iW) {
          gt->signal_weights[wIDs[iW]] = event.genReweight.genParam[iW];
        }
      }
    }

    gt->Fill();

  } // entry loop

  if (DEBUG) { PDebug("PandaLeptonicAnalyzer::Run","Done with entry loop"); }

} // Run()
