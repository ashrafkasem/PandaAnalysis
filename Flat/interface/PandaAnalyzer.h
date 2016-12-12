#ifndef PandaAnalyzer_h
#define PandaAnalyzer_h

// STL
#include "vector"
#include "map"
#include <string>
#include <cmath>

// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLorentzVector.h>

#include "AnalyzerUtilities.h"
#include "GeneralTree.h"

// btag
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
//#include "BTagCalibrationStandalone.h"

// JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

/////////////////////////////////////////////////////////////////////////////
// some misc definitions

/////////////////////////////////////////////////////////////////////////////
// PandaAnalyzer definition
class PandaAnalyzer {
public :
  // configuration enums
  enum PreselectionBit {
   kMonotop     =(1<<0),
   kMonohiggs   =(1<<2),
   kMonojet     =(1<<3),
   kTriggers    =(1<<4),
   kVBF         =(1<<5),
   kRecoil      =(1<<6)
  };

  enum ProcessType { 
    kNone,
    kZ,
    kW,
    kA,
    kTT,
    kTop, // used for non-ttbar top
    kV, // used for non V+jets W or Z
    kH
  };

  enum TriggerBits {
    kMETTrig       =(1<<0),
    kSingleEleTrig =(1<<1),
    kSinglePhoTrig =(1<<2),
    kSingleMuTrig  =(1<<3)
  };

  //////////////////////////////////////////////////////////////////////////////////////

  PandaAnalyzer();
  ~PandaAnalyzer();
  void Init(TTree *tree, TTree *infotree);
  void SetOutputFile(TString fOutName);
  float getMSDcorr(Float_t puppipt, Float_t puppieta);
  void ResetBranches();
  void Run();
  void Terminate();
  bool PassPreselection();
  bool PassEventFilters();
  bool PassGoodLumis();
  void SetDataDir(const char *s);
  void SetPreselectionBit(PreselectionBit b,bool on=true) {
    if (on) 
      preselBits |= b;
    else 
      preselBits &= ~b;
  }

  // public configuration
  void SetFlag(TString flag, bool b=true) { flags[flag]=b; }
  bool isData=false;                         // to do gen matching, etc
  int firstEvent=-1;
  int lastEvent=-1;                          // max events to process; -1=>all
  ProcessType processType=kNone;             // determine what to do the jet matching to

private:

  std::map<TString,bool> flags;

  std::map<panda::PGenParticle*,float> genObjects;         //!< particles we want to match the jets to, and the 'size' of the daughters
  panda::PGenParticle *MatchToGen(double eta, double phi, double r2, int pdgid=0);    //!< private function to match a jet; returns NULL if not found
  std::vector<panda::PObject*> matchPhos, matchEles, matchLeps;
  
  // CMSSW-provided utilities
  FactorizedJetCorrector *ak8MCCorrector=0, *ak8DataCorrector=0; 
//  JetCorrectorParameters *ak8jec=0;
//  JetCorrectionUncertainty *ak8unc=0;

  void calcBJetSFs(TString readername, int flavor, double eta, double pt, double eff, double uncFactor, double &sf, double &sfUp, double &sfDown);
  BTagCalibration *btagCalib=0;
  BTagCalibration *btagCalib_alt=0;
  BTagCalibration *sj_btagCalib;

  std::map<TString,BTagCalibrationReader*> btagReaders; //!< maps "JETTYPE_WP" to a reader 
                                                        // I think we can load multiple flavors in a single reader 
  std::vector<TFile*> openFiles; //!< anything that should be closed
  std::vector<void*> gc; //!< used for misc garbage collection

  // files and histograms containing weights
  TFile *fLepSF=0, *fLepRecoSF=0;
  THCorr2 *hEleVetoLoPU, *hEleTightLoPU;
  THCorr2 *hEleVetoHiPU, *hEleTightHiPU;
  THCorr2 *hMuLooseLoPU, *hMuTightLoPU;
  THCorr2 *hMuLooseHiPU, *hMuTightHiPU;
  THCorr2 *hRecoEleLoPU, *hRecoEleHiPU;
  THCorr2 *hRecoMuLoPU,  *hRecoMuHiPU;
  TFile *fPhoSF=0;
  THCorr2 *hPho=0;

  TFile *fPU=0;
  THCorr1 *hPUWeight;

  TFile *fKFactor=0;
  THCorr1 *hZNLO, *hANLO, *hWNLO;
  THCorr1 *hZLO,  *hALO,  *hWLO;
  THCorr1 *hZEWK, *hAEWK, *hWEWK;
  TFile *fEleTrigB, *fEleTrigE, *fPhoTrig, *fEleTrigLow, *fMetTrig;
  THCorr1 *hEleTrigB, *hEleTrigE, *hPhoTrig, *hMetTrig;
  TFile *fCSVLF, *fCSVHF;
  THCorr1 *hCSVLF, *hCSVHF;
  //THCorr *hEleTrigBUp=0, *hEleTrigBDown=0, *hEleTrigEUp=0, *hEleTrigEDown=0;
  THCorr2 *hEleTrigLow;

  TFile *MSDcorr;
  TF1* puppisd_corrGEN;
  TF1* puppisd_corrRECO_cen;
  TF1* puppisd_corrRECO_for;

  // IO for the analyzer
  TFile *fOut;   // output file is owned by PandaAnalyzer
  TTree *tOut;
  GeneralTree *gt; // essentially a wrapper around tOut
  TTree *tIn=0;  // input tree to read
  unsigned int preselBits=0;

  // objects to read from the tree
  panda::PEvent *event=0;                        // event object
  panda::VGenParticle *genparts=0;               // gen particle objects
  panda::VFatJet *fatjets=0;                     // CA15 fat jets
  panda::VJet *jets=0;                           // AK4 narrow jets
  panda::VElectron *electrons=0;
  panda::VMuon *muons=0;
  panda::VTau *taus=0;
  panda::VPhoton *photons=0;
  panda::PMET *pfmet=0, *puppimet=0;

  // configuration read from output tree
  std::vector<double> betas;
  std::vector<int> Ns; 
  std::vector<int> orders;

};


#endif

