#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TSystem.h"
#include "TMath.h"
#include <algorithm>
#include <vector>
#include "PandaAnalysis/Utilities/src/RoccoR.cc"
#include "PandaAnalysis/Utilities/src/CSVHelper.cc"

#define EGMSCALE 1

using namespace panda;
using namespace std;

PandaAnalyzer::PandaAnalyzer(int debug_/*=0*/) 
{
  DEBUG = debug_;

  if (DEBUG) PDebug("PandaAnalyzer::PandaAnalyzer","Calling constructor");
  gt = new GeneralTree();
  if (DEBUG) PDebug("PandaAnalyzer::PandaAnalyzer","Built GeneralTree");
  ibetas = gt->get_ibetas();
  Ns = gt->get_Ns();
  orders = gt->get_orders();
  if (DEBUG) PDebug("PandaAnalyzer::PandaAnalyzer","Called constructor");
}


PandaAnalyzer::~PandaAnalyzer() 
{
  if (DEBUG) PDebug("PandaAnalyzer::~PandaAnalyzer","Calling destructor");
}


void PandaAnalyzer::ResetBranches() 
{
  genObjects.clear();
  matchPhos.clear();
  matchEles.clear();
  matchLeps.clear();
  looseLeps.clear();
  tightLeps.clear();
  loosePhos.clear();
  cleanedJets.clear();
  isoJets.clear();
  btaggedJets.clear();
  centralJets.clear();
  btagindices.clear();
  genJetsNu.clear();
  fj1 = 0;
  for (TLorentzVector v_ : {vPFMET, vPuppiMET, vpfUW, vpfUZ, vpfUA, vpfU,
                            vpuppiUW, vpuppiUZ, vpuppiUA, vpuppiU,
                            vJet, vBarrelJets})
  {
    v_.SetPtEtaPhiM(0,0,0,0);
  }
  vMETNoMu.SetMagPhi(0,0);
  gt->Reset();
  if (DEBUG) PDebug("PandaAnalyzer::ResetBranches","Reset");
}


void PandaAnalyzer::SetOutputFile(TString fOutName) 
{
  fOut = new TFile(fOutName,"RECREATE");
  fOut->cd();
  tOut = new TTree("events","events");

  fOut->WriteTObject(hDTotalMCWeight);    

  gt->monohiggs      = analysis->monoh;
  gt->vbf            = analysis->vbf;
  gt->fatjet         = analysis->fatjet;
  gt->leptonic       = analysis->complicatedLeptons;
  gt->genPartonStudy = analysis->genPartonStudy;
  gt->btagWeights    = !analysis->btagSFs && analysis->btagWeights;
  gt->useCMVA        = analysis->useCMVA;

  // fill the signal weights
  for (auto& id : wIDs) 
    gt->signal_weights[id] = 1;

  // Build the input tree here 
  gt->WriteTree(tOut);

  if (DEBUG) PDebug("PandaAnalyzer::SetOutputFile","Created output in "+fOutName);
}


int PandaAnalyzer::Init(TTree *t, TH1D *hweights, TTree *weightNames)
{
  if (DEBUG) PDebug("PandaAnalyzer::Init","Starting initialization");
  if (!t || !hweights) {
    PError("PandaAnalyzer::Init","Malformed input!");
    return 0;
  }
  tIn = t;

  event.setStatus(*t, {"!*"}); // turn everything off first

  TString jetname = (analysis->puppi_jets) ? "puppi" : "chs";
  panda::utils::BranchList readlist({"runNumber", "lumiNumber", "eventNumber", "rho", 
                                     "isData", "npv", "npvTrue", "weight", "chsAK4Jets", 
                                     "electrons", "muons", "taus", "photons", 
                                     "pfMet", "caloMet", "puppiMet", "rawMet", 
                                     "recoil","metFilters"});
  readlist.setVerbosity(0);

  if (analysis->ak8)
    readlist += {jetname+"AK8Jets", "subjets", jetname+"AK8Subjets","Subjets"};
  else if (analysis->fatjet) 
    readlist += {jetname+"CA15Jets", "subjets", jetname+"CA15Subjets","Subjets"};
  
  if (analysis->recluster || analysis->bjetRegression)
    readlist.push_back("pfCandidates");

  if (isData) {
    readlist.push_back("triggers");
  } else {
    readlist.push_back("genParticles");
    readlist.push_back("genReweight");
    readlist.push_back("ak4GenJets");
    readlist.push_back("genMet");
  }


  event.setAddress(*t, readlist); // pass the readlist so only the relevant branches are turned on
  if (DEBUG) PDebug("PandaAnalyzer::Init","Set addresses");

  hDTotalMCWeight = new TH1F("hDTotalMCWeight","hDTotalMCWeight",1,0,2);
  hDTotalMCWeight->SetBinContent(1,hweights->GetBinContent(1));

  if (weightNames) {
    if (weightNames->GetEntries()!=377 && weightNames->GetEntries()!=22) {
      PError("PandaAnalyzer::Init",
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
    PError("PandaAnalyzer::Init","This is a signal file, but the weights are missing!");
    return 2;
  }


  // manipulate the output tree
  if (isData) {
    std::vector<TString> droppable = {"mcWeight","scale","scaleUp",
                                      "trueGenBosonPt",
                                      "scaleDown","pdf.*","gen.*","sf_.*"};
    gt->RemoveBranches(droppable,{"sf_phoPurity"});
  }
  if (analysis->genOnly) {
    std::vector<TString> keepable = {"mcWeight","scale","scaleUp",
                                     "scaleDown","pdf.*","gen.*","fj1.*",
                                     "nFatjet","sf_tt.*","sf_qcdTT.*",
                                     "trueGenBosonPt","sf_qcd.*","sf_ewk.*"};
    gt->RemoveBranches({".*"},keepable);
  }

  gt->RemoveBranches({"ak81.*"}); // unused
  
  if (analysis->recluster || analysis->reclusterGen) {
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    activeArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
    areaDef = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*activeArea);
  }

  if (!analysis->fatjet && !analysis->ak8) {
    gt->RemoveBranches({"fj1.*"});
  } else if (analysis->recluster) {
    double radius = 1.5;
    double sdZcut = 0.15;
    double sdBeta = 1.;
    jetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm,radius);
    softDrop = new fastjet::contrib::SoftDrop(sdBeta,sdZcut,radius);
  } else { 
    std::vector<TString> droppable = {"fj1NConst","fj1NSDConst","fj1EFrac100","fj1SDEFrac100"};
    gt->RemoveBranches(droppable);
  }

  if (analysis->reclusterGen) {
    double radius = 0.4;
    jetDefGen = new fastjet::JetDefinition(fastjet::antikt_algorithm,radius);
  }


  if (DEBUG) PDebug("PandaAnalyzer::Init","Finished configuration");

  return 0;
}


panda::GenParticle const *PandaAnalyzer::MatchToGen(double eta, double phi, double radius, int pdgid) 
{
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


void PandaAnalyzer::Terminate() 
{
  fOut->WriteTObject(tOut);
  fOut->Close();

  for (auto *f : fCorrs)
    if (f)
      f->Close();
//   for (unsigned i = 0; i != cN; ++i) {
//     delete h1Corrs[i];
//     printf("1 %i \n", i);
//   }
//   for (unsigned i = 0; i != cN; ++i) {
//     delete h2Corrs[i];
//     printf("2 %i \n", i);
//   }
  for (auto *h : h1Corrs)
    delete h;
  for (auto *h : h2Corrs)
    delete h;

  delete btagCalib;
  delete sj_btagCalib;
  for (auto *reader : btagReaders )
    delete reader;

  for (auto& iter : ak8UncReader)
    delete iter.second;

  delete ak8JERReader;

  for (auto& iter : ak4UncReader)
    delete iter.second;

  for (auto& iter : ak4ScaleReader) {
    printf("trying to delete: |%s| at %p\n",iter.first.Data(),iter.second);
    delete iter.second;
  }

  delete ak4JERReader;

  delete activeArea;
  delete areaDef;
  delete jetDef;
  delete jetDefGen;
  delete softDrop;

  delete hDTotalMCWeight;
  if (DEBUG) PDebug("PandaAnalyzer::Terminate","Finished with output");
}

void PandaAnalyzer::OpenCorrection(CorrectionType ct, TString fpath, TString hname, int dim) 
{
  fCorrs[ct] = TFile::Open(fpath);
  if (dim==1) 
    h1Corrs[ct] = new THCorr1((TH1D*)fCorrs[ct]->Get(hname));
  else if (dim==2)
    h2Corrs[ct] = new THCorr2((TH2D*)fCorrs[ct]->Get(hname));
  else 
    f1Corrs[ct] = new TF1Corr((TF1*)fCorrs[ct]->Get(hname));
}

double PandaAnalyzer::GetCorr(CorrectionType ct, double x, double y) 
{
  if (h1Corrs[ct]!=0) {
    return h1Corrs[ct]->Eval(x); 
  } else if (h2Corrs[ct]!=0) {
    return h2Corrs[ct]->Eval(x,y);
  } else if (f1Corrs[ct]!=0) {
    return f1Corrs[ct]->Eval(x);
  } else {
    PError("PandaAnalyzer::GetCorr",
       TString::Format("No correction is defined for CorrectionType=%u",ct));
    return 1;
  }
}
double PandaAnalyzer::GetError(CorrectionType ct, double x, double y) {
  if (h1Corrs[ct]!=0) {
    return h1Corrs[ct]->Error(x); 
  } else if (h2Corrs[ct]!=0) {
    return h2Corrs[ct]->Error(x,y);
  } else {
    PError("PandaAnalyzer::GetError",
       TString::Format("No correction is defined for CorrectionType=%u",ct));
    return 1;
  }
}
void PandaAnalyzer::SetDataDir(const char *s) 
{
  TString dirPath(s);
  dirPath += "/";

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Starting loading of data");

  // pileup
  OpenCorrection(cNPV,dirPath+"moriond17/normalized_npv.root","data_npv_Wmn",1);
  OpenCorrection(cPU,dirPath+"moriond17/puWeights_80x_37ifb.root","puWeights",1);
  OpenCorrection(cPUUp,dirPath+"moriond17/puWeights_80x_37ifb.root","puWeightsUp",1);
  OpenCorrection(cPUDown,dirPath+"moriond17/puWeights_80x_37ifb.root","puWeightsDown",1);

  if (analysis->complicatedLeptons) {
    // Corrections checked out from Gui's repository on Nov 12, 2017 ~DGH
    // https://github.com/GuillelmoGomezCeballos/MitAnalysisRunII/tree/master/data/80x
    OpenCorrection(cZHEwkCorr,dirPath+"leptonic/Zll_nloEWK_weight_unnormalized.root","SignalWeight_nloEWK_rebin",1);
    OpenCorrection(cZHEwkCorrUp  ,dirPath+"leptonic/Zll_nloEWK_weight_unnormalized.root","SignalWeight_nloEWK_up_rebin",1);
    OpenCorrection(cZHEwkCorrDown,dirPath+"leptonic/Zll_nloEWK_weight_unnormalized.root","SignalWeight_nloEWK_down_rebin",1);
    OpenCorrection(cMuLooseID,dirPath+"leptonic/muon_scalefactors_37ifb.root","scalefactors_MuonLooseId_Muon",2);
    OpenCorrection(cMuMediumID,dirPath+"leptonic/scalefactors_80x_dylan_37ifb.root","scalefactors_Medium_Muon",2);
    OpenCorrection(cMuTightID,dirPath+"leptonic/muon_scalefactors_37ifb.root","scalefactors_TightId_Muon",2);
    OpenCorrection(cMuLooseIso,dirPath+"leptonic/muon_scalefactors_37ifb.root","scalefactors_Iso_MuonLooseId",2);
    OpenCorrection(cMuMediumIso,dirPath+"leptonic/muon_scalefactors_37ifb.root","scalefactors_Iso_MuonMediumId",2);
    OpenCorrection(cMuTightIso,dirPath+"leptonic/muon_scalefactors_37ifb.root","scalefactors_Iso_MuonTightId",2);
    OpenCorrection(cMuReco,dirPath+"leptonic/Tracking_EfficienciesAndSF_BCDEFGH.root","ratio_eff_eta3_dr030e030_corr",1);
    OpenCorrection(cEleVeto,dirPath+"moriond17/scaleFactor_electron_summer16.root","scaleFactor_electron_vetoid_RooCMSShape_pu_0_100",2);
    OpenCorrection(cEleLoose,dirPath+"leptonic/scalefactors_80x_egpog_37ifb.root","scalefactors_Loose_Electron",2);
    OpenCorrection(cEleMedium,dirPath+"leptonic/scalefactors_80x_dylan_37ifb.root","scalefactors_Medium_Electron",2);
    OpenCorrection(cEleTight,dirPath+"leptonic/scalefactors_80x_egpog_37ifb.root","scalefactors_Tight_Electron",2);
    OpenCorrection(cEleReco,dirPath+"leptonic/scalefactors_80x_egpog_37ifb.root","scalefactors_Reco_Electron",2);
    // EWK corrections 
    OpenCorrection(cWZEwkCorr,dirPath+"leptonic/data.root","hEWKWZCorr",1);
    OpenCorrection(cqqZZQcdCorr,dirPath+"leptonic/data.root","hqqZZKfactor",2);
  } else {
    OpenCorrection(cEleVeto,dirPath+"moriond17/scaleFactor_electron_summer16.root","scaleFactor_electron_vetoid_RooCMSShape_pu_0_100",2);
    OpenCorrection(cEleTight,dirPath+"moriond17/scaleFactor_electron_summer16.root","scaleFactor_electron_tightid_RooCMSShape_pu_0_100",2);
    OpenCorrection(cEleReco,dirPath+"moriond17/scaleFactor_electron_reco_summer16.root","scaleFactor_electron_reco_RooCMSShape_pu_0_100",2);
    OpenCorrection(cMuLooseID,dirPath+"moriond17/muon_scalefactors_37ifb.root","scalefactors_MuonLooseId_Muon",2);
    OpenCorrection(cMuLooseIso,dirPath+"moriond17/muon_scalefactors_37ifb.root","scalefactors_Iso_MuonLooseId",2);
    OpenCorrection(cMuTightID,dirPath+"moriond17/muon_scalefactors_37ifb.root","scalefactors_TightId_Muon",2);
    OpenCorrection(cMuTightIso,dirPath+"moriond17/muon_scalefactors_37ifb.root","scalefactors_Iso_MuonTightId",2);
    OpenCorrection(cMuReco,dirPath+"moriond17/Tracking_12p9.root","htrack2",1);
  }
  // photons
  OpenCorrection(cPho,dirPath+"moriond17/scalefactors_80x_medium_photon_37ifb.root",
                 "EGamma_SF2D",2);

  // triggers
  OpenCorrection(cTrigMET,dirPath+"moriond17/metTriggerEfficiency_recoil_monojet_TH1F.root",
                 "hden_monojet_recoil_clone_passed",1);
  OpenCorrection(cTrigEle,dirPath+"moriond17/eleTrig.root","hEffEtaPt",2);
  OpenCorrection(cTrigPho,dirPath+"moriond17/photonTriggerEfficiency_photon_TH1F.root",
                 "hden_photonpt_clone_passed",1);
  OpenCorrection(cTrigMETZmm,dirPath+"moriond17/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root",
                 "hden_monojet_recoil_clone_passed",1);

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded scale factors");

  // kfactors
  TFile *fKFactor = 0;
  if (analysis->vbf)
    fKFactor = new TFile(dirPath+"vbf16/kqcd/kfactor_24bins.root"); 
  else
    fKFactor = new TFile(dirPath+"kfactors.root"); 
  fCorrs[cZNLO] = fKFactor; // just for garbage collection

  TH1D *hZLO    = (TH1D*)fKFactor->Get("ZJets_LO/inv_pt");
  TH1D *hWLO    = (TH1D*)fKFactor->Get("WJets_LO/inv_pt");
  TH1D *hALO    = (TH1D*)fKFactor->Get("GJets_LO/inv_pt_G");

  h1Corrs[cZNLO] = new THCorr1((TH1D*)fKFactor->Get("ZJets_012j_NLO/nominal"));
  h1Corrs[cWNLO] = new THCorr1((TH1D*)fKFactor->Get("WJets_012j_NLO/nominal"));
  h1Corrs[cANLO] = new THCorr1((TH1D*)fKFactor->Get("GJets_1j_NLO/nominal_G"));

  h1Corrs[cZEWK] = new THCorr1((TH1D*)fKFactor->Get("EWKcorr/Z"));
  h1Corrs[cWEWK] = new THCorr1((TH1D*)fKFactor->Get("EWKcorr/W"));
  h1Corrs[cAEWK] = new THCorr1((TH1D*)fKFactor->Get("EWKcorr/photon"));

  h1Corrs[cZEWK]->GetHist()->Divide(h1Corrs[cZNLO]->GetHist());     
  h1Corrs[cWEWK]->GetHist()->Divide(h1Corrs[cWNLO]->GetHist());     
  h1Corrs[cAEWK]->GetHist()->Divide(h1Corrs[cANLO]->GetHist());

  h1Corrs[cZNLO]->GetHist()->Divide(hZLO);    
  h1Corrs[cWNLO]->GetHist()->Divide(hWLO);    
  h1Corrs[cANLO]->GetHist()->Divide(hALO);

  OpenCorrection(cANLO2j,dirPath+"moriond17/histo_photons_2jet.root","Func",1);

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded k factors");

  if (analysis->vbf) {

    OpenCorrection(cVBF_ZNLO,dirPath+"vbf16/kqcd/mjj/merged_zvv.root","h_kfactors_shape",2);
    OpenCorrection(cVBF_WNLO,dirPath+"vbf16/kqcd/mjj/merged_wlv.root","h_kfactors_shape",2);
    OpenCorrection(cVBF_ZllNLO,dirPath+"vbf16/kqcd/mjj/merged_zll.root","h_kfactors_shape",2);

    OpenCorrection(cVBFTight_ZNLO,dirPath+"vbf16/kqcd/mjj/merged_zvv.root","h_kfactors_cc",1);
    OpenCorrection(cVBFTight_WNLO,dirPath+"vbf16/kqcd/mjj/merged_wlv.root","h_kfactors_cc",1);
    OpenCorrection(cVBFTight_ZllNLO,dirPath+"vbf16/kqcd/mjj/merged_zll.root","h_kfactors_cc",1);

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded VBF k factors");
    /*
    TFile *fKFactor_VBFZ = new TFile(dirPath+"vbf16/kqcd/kfactor_VBF_zjets_v2.root");
    h1Corrs[cVBF_ZNLO] = new THCorr1((TH1D*)fKFactor_VBFZ->Get("bosonPt_NLO_vbf_relaxed"));
    h1Corrs[cVBF_ZNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFZ->Get("bosonPt_LO_vbf_relaxed"));
    h1Corrs[cVBF_ZNLO]->GetHist()->Multiply((TH1D*)fKFactor_VBFZ->Get("bosonPt_LO_monojet"));
    h1Corrs[cVBF_ZNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFZ->Get("bosonPt_NLO_monojet"));
    h1Corrs[cVBFTight_ZNLO] = new THCorr1((TH1D*)fKFactor_VBFZ->Get("bosonPt_NLO_vbf"));
    h1Corrs[cVBFTight_ZNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFZ->Get("bosonPt_LO_vbf"));
    h1Corrs[cVBFTight_ZNLO]->GetHist()->Multiply((TH1D*)fKFactor_VBFZ->Get("bosonPt_LO_monojet"));
    h1Corrs[cVBFTight_ZNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFZ->Get("bosonPt_NLO_monojet"));

    TFile *fKFactor_VBFW = new TFile(dirPath+"vbf16/kqcd/kfactor_VBF_wjets_v2.root");
    h1Corrs[cVBF_WNLO] = new THCorr1((TH1D*)fKFactor_VBFW->Get("bosonPt_NLO_vbf_relaxed"));
    h1Corrs[cVBF_WNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFW->Get("bosonPt_LO_vbf_relaxed"));
    h1Corrs[cVBF_WNLO]->GetHist()->Multiply((TH1D*)fKFactor_VBFW->Get("bosonPt_LO_monojet"));
    h1Corrs[cVBF_WNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFW->Get("bosonPt_NLO_monojet"));
    h1Corrs[cVBFTight_WNLO] = new THCorr1((TH1D*)fKFactor_VBFW->Get("bosonPt_NLO_vbf"));
    h1Corrs[cVBFTight_WNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFW->Get("bosonPt_LO_vbf"));
    h1Corrs[cVBFTight_WNLO]->GetHist()->Multiply((TH1D*)fKFactor_VBFW->Get("bosonPt_LO_monojet"));
    h1Corrs[cVBFTight_WNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFW->Get("bosonPt_NLO_monojet"));
    */

    OpenCorrection(cVBF_EWKZ,dirPath+"vbf16/kewk/kFactor_ZToNuNu_pT_Mjj.root",
                   "TH2F_kFactor",2);
    OpenCorrection(cVBF_EWKW,dirPath+"vbf16/kewk/kFactor_WToLNu_pT_Mjj.root",
                   "TH2F_kFactor",2);

    OpenCorrection(cVBF_TrigMET,dirPath+"vbf16/trig/fit_nmu1.root",
                   "f_eff",3);
    OpenCorrection(cVBF_TrigMETZmm,dirPath+"vbf16/trig/fit_nmu2.root",
                   "f_eff",3);

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded VBF k factors");
  }

  OpenCorrection(cBadECALJets,dirPath+"vbf16/hotjets-runBCDEFGH.root",
                 "h2jet",2);


  if (analysis->btagSFs) {
    // btag SFs
    btagCalib = new BTagCalibration("csvv2",(dirPath+"moriond17/CSVv2_Moriond17_B_H.csv").Data());
    btagReaders[bJetL] = new BTagCalibrationReader(BTagEntry::OP_LOOSE,"central",{"up","down"});
    btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_B,"comb");
    btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_C,"comb");
    btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_UDSG,"incl");

    sj_btagCalib = new BTagCalibration("csvv2",(dirPath+"moriond17/subjet_CSVv2_Moriond17_B_H.csv").Data());
    btagReaders[bSubJetL] = new BTagCalibrationReader(BTagEntry::OP_LOOSE,"central",{"up","down"});
    btagReaders[bSubJetL]->load(*sj_btagCalib,BTagEntry::FLAV_B,"lt");
    btagReaders[bSubJetL]->load(*sj_btagCalib,BTagEntry::FLAV_C,"lt");
    btagReaders[bSubJetL]->load(*sj_btagCalib,BTagEntry::FLAV_UDSG,"incl");

    btagReaders[bJetM] = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"});
    btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_B,"comb");
    btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_C,"comb");
    btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_UDSG,"incl");

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded btag SFs");
  } else if(analysis->btagWeights) {
    if (analysis->useCMVA) 
      cmvaReweighter = new CSVHelper("PandaAnalysis/data/csvweights/cmva_rwt_fit_hf_v0_final_2017_3_29.root"   , "PandaAnalysis/data/csvweights/cmva_rwt_fit_lf_v0_final_2017_3_29.root"   , 5);
    else
      csvReweighter  = new CSVHelper("PandaAnalysis/data/csvweights/csv_rwt_fit_hf_v2_final_2017_3_29test.root", "PandaAnalysis/data/csvweights/csv_rwt_fit_lf_v2_final_2017_3_29test.root", 5);
  }

  // bjet regression
  if (analysis->bjetRegression){
    bjetreg_vars = new float[10];

    bjetreg_reader->AddVariable("jetPt[hbbjtidx[0]]",&bjetreg_vars[0]);
    bjetreg_reader->AddVariable("nJot",&bjetreg_vars[1]);
    bjetreg_reader->AddVariable("jetEta[hbbjtidx[0]]",&bjetreg_vars[2]);
    bjetreg_reader->AddVariable("jetE[hbbjtidx[0]]",&bjetreg_vars[3]);
    bjetreg_reader->AddVariable("npv",&bjetreg_vars[4]);
    bjetreg_reader->AddVariable("jetLeadingTrkPt[hbbjtidx[0]]",&bjetreg_vars[5]);
    bjetreg_reader->AddVariable("jetLeadingLepPt[hbbjtidx[0]]",&bjetreg_vars[6]);
    bjetreg_reader->AddVariable("jetNLep[hbbjtidx[0]]",&bjetreg_vars[7]);
    bjetreg_reader->AddVariable("jetEMFrac[hbbjtidx[0]]",&bjetreg_vars[8]);
    bjetreg_reader->AddVariable("jetHadFrac[hbbjtidx[0]]",&bjetreg_vars[9]);

    bjetreg_reader->BookMVA( "BDT method", dirPath+"trainings/bjet_regression_v0.weights.xml" );    

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded bjet regression weights");
  }


  if (analysis->monoh) {
    // mSD corr
    MSDcorr = new TFile(dirPath+"/puppiCorr.root");
    puppisd_corrGEN = (TF1*)MSDcorr->Get("puppiJECcorr_gen");;
    puppisd_corrRECO_cen = (TF1*)MSDcorr->Get("puppiJECcorr_reco_0eta1v3");
    puppisd_corrRECO_for = (TF1*)MSDcorr->Get("puppiJECcorr_reco_1v3eta2v5");

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded mSD correction");
  }

  if (analysis->rerunJES) {
    TString jecV = "V4", jecReco = "23Sep2016"; 
    TString jecVFull = jecReco+jecV;
    ak8UncReader["MC"] = new JetCorrectionUncertainty(
       (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_Uncertainty_AK8PFPuppi.txt").Data()
      );
    std::vector<TString> eraGroups = {"BCD","EF","G","H"};
    for (auto e : eraGroups) {
      ak8UncReader["data"+e] = new JetCorrectionUncertainty(
         (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_Uncertainty_AK8PFPuppi.txt").Data()
        );
    }

    ak8JERReader = new JERReader(dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_SF_AK8PFPuppi.txt",
                                 dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt");


    ak4UncReader["MC"] = new JetCorrectionUncertainty(
       (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_Uncertainty_AK4PFPuppi.txt").Data()
      );
    for (auto e : eraGroups) {
      ak4UncReader["data"+e] = new JetCorrectionUncertainty(
         (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_Uncertainty_AK4PFPuppi.txt").Data()
        );
    }

    ak4JERReader = new JERReader(dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_SF_AK4PFPuppi.txt",
                                 dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_PtResolution_AK4PFPuppi.txt");

    std::vector<JetCorrectorParameters> params = {
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L1FastJet_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L2Relative_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L3Absolute_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L2L3Residual_AK4PFPuppi.txt").Data())
    };
    ak4ScaleReader["MC"] = new FactorizedJetCorrector(params);
    if (DEBUG>1) PDebug("PandaAnalyzer::SetDataDir","Loaded JES for AK4 MC");
    for (auto e : eraGroups) {
      params = {
        JetCorrectorParameters(
          (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L1FastJet_AK4PFPuppi.txt").Data()),
        JetCorrectorParameters(
          (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L2Relative_AK4PFPuppi.txt").Data()),
        JetCorrectorParameters(
          (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L3Absolute_AK4PFPuppi.txt").Data()),
        JetCorrectorParameters(
          (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L2L3Residual_AK4PFPuppi.txt").Data())
      };
      ak4ScaleReader["data"+e] = new FactorizedJetCorrector(params);
      if (DEBUG>1) PDebug("PandaAnalyzer::SetDataDir","Loaded JES for AK4 "+e);
    }

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded JES/R");
  }

}


void PandaAnalyzer::AddGoodLumiRange(int run, int l0, int l1) 
{
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) { // don't know about this run yet
    std::vector<LumiRange> newLumiList;
    newLumiList.emplace_back(l0,l1);
    goodLumis[run] = newLumiList;
  } else {
    run_->second.emplace_back(l0,l1);
  }
}


bool PandaAnalyzer::PassGoodLumis(int run, int lumi) 
{
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) {
    // matched no run
    if (DEBUG) 
      PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i",run));
    return false;
  }

  // found the run, now look for a lumi range
  for (auto &range : run_->second) {
    if (range.Contains(lumi)) {
      if (DEBUG) 
        PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Accepting run=%i, lumi=%i",run,lumi));
      return true;
    }
  }

  // matched no lumi range
  if (DEBUG) 
    PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i, lumi=%i",run,lumi));
  return false;
}


bool PandaAnalyzer::PassPreselection() 
{
  // TODO: refactor this function
  // was originally written this way to handle more complex conditions
  // like triggers, but could probably clean it up with a Condition class
  
  if (preselBits==0)
    return true;
  bool isGood=false;

  if (preselBits & kGenBosonPt) {
    if (gt->trueGenBosonPt > 100)
      isGood = true; 
  }

  if (preselBits & kFatjet) {
    if (gt->fj1Pt>250)
      isGood = true;
  }

  float max_puppi = std::max({gt->puppimet, gt->puppiUZmag, gt->puppiUWmag, gt->puppiUAmag});
  float max_pf = std::max({gt->pfmet, gt->pfUZmag, gt->pfUWmag, gt->pfUAmag});
  float max_pfUp = std::max({gt->pfmetUp, gt->pfUZmagUp, gt->pfUWmagUp, gt->pfUAmagUp});
  float max_pfDown = std::max({gt->pfmetDown, gt->pfUZmagDown, gt->pfUWmagDown, gt->pfUAmagDown});

  if (preselBits & kRecoil) {
    if ( max_pfDown>200 || max_pf>200 || max_pfUp>200 || max_puppi>200 ) {
      isGood = true;
    }
  }
  if (preselBits & kRecoil50) {
    if ( gt->pfmet>100 ) { // this will never cause any confusion, I'm sure
      isGood = true;
    }
  }
  if (preselBits & kMonotop) {
    if (gt->nFatjet>=1 && gt->fj1Pt>200) {
      if ( max_pf>200 || max_puppi>200) {
        isGood = true;
      }
    }
  }
  if (preselBits & kMonojet) {
    if (true) {
      if ( max_pfDown>200 || max_pf>200 || max_pfUp>200 || max_puppi>200 ) {
        isGood = true;
      }
    }
  }
  if (preselBits & kMonohiggs) {
    if ((gt->nFatjet>=1 && gt->fj1Pt>200) || gt->hbbpt>150 ) {
      if ( max_pf>175 || max_puppi>175) {
        isGood = true;
      }
    }
  }

  if(preselBits & kVHBB) {
    // ZnnHbb
    if (
      gt->pfmet>150 && 
      gt->nJet>=2 && gt->jetPt[0]>50 && gt->jetPt[1]>50 &&
      (gt->hbbpt>50 || (gt->nFatjet>0 && gt->fj1Pt>200))
    ) isGood=true;
    // WlnHbb
    else if (
      gt->pfmet>30 &&
      gt->nJet>=2 && gt->jetPt[0]>25 && gt->jetPt[1]>25 &&
      (
       (gt->nTightElectron >0 && gt->electronPt[0]>25) ||
       (gt->nTightMuon > 0 && gt->muonPt[0]>25)
      ) &&
      (gt->hbbpt>50 || (gt->nFatjet>0 && gt->fj1Pt>200))
    ) isGood=true;
    // ZllHbb
    else if (
      gt->nJet>=2 && gt->jetPt[0]>25 && gt->jetPt[1]>25 &&
      (
       (
        gt->nTightElectron>0 && 
        gt->nLooseElectron>1 &&
        gt->electronPt[0]>25 && 
        gt->electronPt[1]>20
       ) || (
        gt->nTightMuon > 0 && 
        gt->nLooseMuon>1 &&
        gt->muonPt[0]>25 &&
        gt->muonPt[1]>20 
       )
      ) &&
      (gt->hbbpt>50 || (gt->nFatjet>0 && gt->fj1Pt>200))
    ) isGood=true;
  }
  // anded with the rest
  if (preselBits & kPassTrig) {
    isGood = (!isData) || (gt->trigger != 0);
  }

  tr->TriggerEvent("presel");
  
  return isGood;
}



// run
void PandaAnalyzer::Run() 
{

  fOut->cd(); // to be absolutely sure

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
  genBosonPtMin=150, genBosonPtMax=1000;
  if (!isData && h1Corrs[cZNLO]) {
    genBosonPtMin = h1Corrs[cZNLO]->GetHist()->GetBinCenter(1);
    genBosonPtMax = h1Corrs[cZNLO]->GetHist()->GetBinCenter(h1Corrs[cZNLO]->GetHist()->GetNbinsX());
  }

  if (analysis->ak8) {
    if (analysis->puppi_jets)
      fatjets = &event.puppiAK8Jets;
    else
      fatjets = &event.chsAK8Jets;
  } else if (analysis->fatjet) {
    if (analysis->puppi_jets)
      fatjets = &event.puppiCA15Jets;
    else
      fatjets = &event.chsCA15Jets;
  }

  jets = &event.chsAK4Jets;

  // these are bins of b-tagging eff in pT and eta, derived in 8024 TT MC
  // TODO: don't hardcode these 
  std::vector<double> vbtagpt {20.0,50.0,80.0,120.0,200.0,300.0,400.0,500.0,700.0,1000.0};
  std::vector<double> vbtageta {0.0,0.5,1.5,2.5};
  lfeff  = {{0.081,0.065,0.060,0.063,0.072,0.085,0.104,0.127,0.162},
            {0.116,0.097,0.092,0.099,0.112,0.138,0.166,0.185,0.222},
            {0.173,0.145,0.149,0.175,0.195,0.225,0.229,0.233,0.250}};
  ceff = {{0.377,0.389,0.391,0.390,0.391,0.375,0.372,0.392,0.435},
          {0.398,0.407,0.416,0.424,0.424,0.428,0.448,0.466,0.500},
          {0.375,0.389,0.400,0.425,0.437,0.459,0.481,0.534,0.488}};
  beff = {{0.791,0.815,0.825,0.835,0.821,0.799,0.784,0.767,0.760},
          {0.794,0.816,0.829,0.836,0.823,0.804,0.798,0.792,0.789},
          {0.739,0.767,0.780,0.789,0.776,0.771,0.779,0.787,0.806}};
  btagpt = Binner(vbtagpt);
  btageta = Binner(vbtageta);
  if (analysis->complicatedLeptons) {
    if (DEBUG) PDebug("PandaAnalyzer::Run","Loading the Rochester corrections with random seed 3393");
    // TO DO: Hard coded to 2016 rochester corrections for now, need to do this in a better way later
    TString dirPath1 = TString(gSystem->Getenv("CMSSW_BASE")) + "/src/";
    rochesterCorrection = new RoccoR(Form("%sPandaAnalysis/data/rcdata.2016.v3",dirPath1.Data()));
    rng=TRandom3(3393); //Dylan's b-day
  }



  std::vector<unsigned int> metTriggers;
  std::vector<unsigned int> eleTriggers;
  std::vector<unsigned int> phoTriggers;
  std::vector<unsigned int> muTriggers;
  std::vector<unsigned int> jetTriggers;

  if (isData) {
    if (DEBUG) PDebug("PandaAnalyzer::Run","Loading the trigger paths");
    std::vector<TString> paths;
    paths = {
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
    triggerHandlers[kMETTrig].addTriggers(paths);

    if (analysis->complicatedLeptons)
      paths = {
            "HLT_Ele27_WPTight_Gsf",
      };
    else
      paths = {
            "HLT_Ele27_WP85_Gsf",
            "HLT_Ele27_WPLoose_Gsf",
            "HLT_Ele105_CaloIdVT_GsfTrkIdT",
            "HLT_Ele27_WPTight_Gsf",
            "HLT_Ele30_WPTight_Gsf",
            "HLT_Ele27_eta2p1_WPTight_Gsf",
            "HLT_Ele32_eta2p1_WPTight_Gsf",
            "HLT_Ele35_WPLoose_Gsf",
            "HLT_ECALHT800"
      };
    
    triggerHandlers[kSingleEleTrig].addTriggers(paths);
    
    paths = {
	      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
	      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
	      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
	      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"
    };
    triggerHandlers[kDoubleMuTrig].addTriggers(paths);
    paths = {
	      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
	      "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf"
    };
    triggerHandlers[kDoubleEleTrig].addTriggers(paths);
    
    paths = {
	      "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
	      "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
	      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
	      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
	      "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
	      "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
	      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
	      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"
    };
    triggerHandlers[kEMuTrig].addTriggers(paths);

    paths = {
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
    triggerHandlers[kSinglePhoTrig].addTriggers(paths);

    paths = {
          "HLT_PFHT650",
          "HLT_PFHT900",
          "HLT_PFJet500",
          "HLT_PFJet450",
          "HLT_PFJet320",
    };
    triggerHandlers[kJetHTTrig].addTriggers(paths);

    if (analysis->complicatedLeptons)
      paths = {
            "HLT_IsoMu24",
            "HLT_IsoTkMu24"
      };
    else
      paths = {
            "HLT_IsoMu20",
            "HLT_IsoMu22",
            "HLT_IsoMu24",
      };
    triggerHandlers[kSingleMuTrig].addTriggers(paths);
    
    RegisterTriggers();
  }

  if (analysis->ak8)
    FATJETMATCHDR2 = 0.64;

  fOut->cd(); // to be absolutely sure

  // set up reporters
  unsigned int iE=0;
  ProgressReporter pr("PandaAnalyzer::Run",&iE,&nEvents,10);
  tr = new TimeReporter("PandaAnalyzer::Run",DEBUG);


  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    tr->Start();
    pr.Report();
    ResetBranches();
    event.getEntry(*tIn,iE);


    tr->TriggerEvent(TString::Format("GetEntry %u",iE));
    if (DEBUG>2) {
      PDebug("PandaAnalyzer::Run::Dump","");
      event.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.photons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.muons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.electrons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.chsAK4Jets.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.pfMet.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.metMuOnlyFix.print(std::cout, 2);
      std::cout << std::endl;
    }

    if (!RecoilPresel())
      continue;

    // event info
    gt->mcWeight = event.weight;
    gt->runNumber = event.runNumber;
    gt->lumiNumber = event.lumiNumber;
    gt->eventNumber = event.eventNumber;
    gt->npv = event.npv;
    gt->pu = event.npvTrue;
    gt->metFilter = (event.metFilters.pass()) ? 1 : 0;
    gt->metFilter = (gt->metFilter==1 && !event.metFilters.badPFMuons) ? 1 : 0;
    gt->metFilter = (gt->metFilter==1 && !event.metFilters.badChargedHadrons) ? 1 : 0;
    if (isData) {
      // check the json
      if (!PassGoodLumis(gt->runNumber,gt->lumiNumber))
        continue;

      // save triggers
      for (unsigned iT = 0; iT != kNTrig; ++iT) {
        auto &th = triggerHandlers.at(iT);
        for (auto iP : th.indices) {
          if (event.triggerFired(iP)) {
              gt->trigger |= (1 << iT);
              break;
          }
        }
      }
    } else { // !isData
      gt->sf_npv = GetCorr(cNPV,gt->npv);
      gt->sf_pu = GetCorr(cPU,gt->pu);
      gt->sf_puUp = GetCorr(cPUUp,gt->pu);
      gt->sf_puDown = GetCorr(cPUDown,gt->pu);
    }

    if (analysis->rerunJES)
      SetupJES();

    tr->TriggerEvent("initialize");

    // met
    gt->pfmetRaw = event.rawMet.pt;
    gt->pfmet = event.pfMet.pt;
    gt->pfmetphi = event.pfMet.phi;
    gt->calomet = event.caloMet.pt;
    gt->sumETRaw = event.pfMet.sumETRaw;
    gt->puppimet = event.puppiMet.pt;
    gt->puppimetphi = event.puppiMet.phi;
    GetMETSignificance();
    vPFMET.SetPtEtaPhiM(gt->pfmet,0,gt->pfmetphi,0);
    vPuppiMET.SetPtEtaPhiM(gt->puppimet,0,gt->puppimetphi,0);
    vMETNoMu.SetMagPhi(gt->pfmet,gt->pfmetphi); //       for trigger eff
    if (analysis->varyJES) {
      gt->pfmetUp = event.pfMet.ptCorrUp;
      gt->pfmetDown = event.pfMet.ptCorrDown;
    }

    tr->TriggerEvent("met");

    // electrons and muons
    if (analysis->complicatedLeptons) {
      ComplicatedLeptons();
      GenStudyEWK();
    } else SimpleLeptons();
    
    // photons
    Photons();

    // recoil!
    if (analysis->recoil)
      Recoil();

    // fatjets
    if (analysis->fatjet) {
      FatjetBasics();
      if (analysis->recluster)
        FatjetRecluster();
      tr->TriggerEvent("fatjet");
    }

    // first identify interesting jets
    JetBasics();

    if (analysis->monoh) {
      // Higgs reconstruction for resolved analysis - highest pt pair of b jets
      JetHbbReco();
    }

    Taus();

    if (!analysis->genOnly && !PassPreselection()) // only check reco presel here
      continue;

    gt->sf_tt = 1; gt->sf_tt_ext = 1; gt->sf_tt_bound = 1;
    gt->sf_tt8TeV = 1; gt->sf_tt8TeV_ext = 1; gt->sf_tt8TeV_bound = 1;
    gt->sf_qcdTT = 1;
    gt->sf_qcdV=1; gt->sf_ewkV=1;
    gt->sf_qcdV_VBF=1;
    gt->sf_qcdV_VBFTight=1;
    gt->sf_qcdV_VBF2l=1;
    gt->sf_qcdV_VBF2lTight=1;
    gt->sf_eleTrig=1; gt->sf_metTrig=1; gt->sf_phoTrig=1; gt->sf_metTrigZmm=1;
    gt->sf_metTrigVBF=1; gt->sf_metTrigZmmVBF=1;
    gt->sf_metTrigVBF=1; gt->sf_metTrigZmmVBF=1;
    gt->sf_lepID=1; gt->sf_lepIso=1; gt->sf_lepTrack=1;
    gt->sf_pho=1;
    gt->scaleUp = 1; gt->scaleDown = 1;
    gt->pdfUp = 1; gt->pdfDown = 1;

    if (!isData) {
      if (analysis->fatjet)
        FatjetMatching();
      else if(analysis->genPartonStudy)
        GenPartonStudy();

      if (analysis->btagSFs)
        JetBtagSFs();
      else if(analysis->btagWeights) {
        JetCMVAWeights();
      }

      TopPTReweight();
      VJetsReweight();

      if (analysis->vbf)
        SaveGenLeptons();

      TriggerEffs();

      SignalInfo();

      if(!analysis->complicatedLeptons) LeptonSFs();

      PhotonSFs();

      QCDUncs();
      SignalReweights();

      if (analysis->reclusterGen && analysis->monoh) {
        GenJetsNu();
        MatchGenJets(genJetsNu);
      }
    }

    
    if (analysis->genOnly && !PassPreselection()) // only check gen presel here
      continue;

    gt->Fill();

  } // entry loop

  if (DEBUG) { PDebug("PandaAnalyzer::Run","Done with entry loop"); }

} // Run()

