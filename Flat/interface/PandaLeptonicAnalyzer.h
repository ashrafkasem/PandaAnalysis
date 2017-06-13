#ifndef PandaLeptonicAnalyzer_h
#define PandaLeptonicAnalyzer_h

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
#include "GeneralLeptonicTree.h"

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
// PandaLeptonicAnalyzer definition
class PandaLeptonicAnalyzer {
public :
    // configuration enums
    enum PreselectionBit {
     kTriggers   =(1<<0),
     kVBF        =(1<<1),
     kRecoil     =(1<<2),
     kLepton     =(1<<3)
    };

    enum SelectionBit {
     kLoose   =(1<<0),
     kFake    =(1<<1),
     kMedium  =(1<<2),
     kTight   =(1<<3)
    };

    enum ProcessType { 
        kNone,
        kZ,
        kW,
        kA,
        kZEWK,
        kWEWK,
        kTT,
        kTop, // used for non-ttbar top
        kV, // used for non V+jets W or Z
        kH,
        kSignal,
    };

    enum TriggerBits {
        kMETTrig       =(1<<0),
        kSinglePhoTrig =(1<<1),
        kMuEGTrig      =(1<<2),
        kMuMuTrig      =(1<<3),
        kMuTrig        =(1<<4),
        kEGEGTrig      =(1<<5),
        kEGTrig        =(1<<6)
    };

    //////////////////////////////////////////////////////////////////////////////////////

    PandaLeptonicAnalyzer(int debug_=0);
    ~PandaLeptonicAnalyzer();
    int Init(TTree *tree, TH1D *hweights, TTree *weightNames=0);
    void SetOutputFile(TString fOutName);
    void ResetBranches();
    void Run();
    void Terminate();
    void SetDataDir(const char *s2);
    void SetPreselectionBit(PreselectionBit b,bool on=true) {
        if (on) 
            preselBits |= b;
        else 
            preselBits &= ~b;
    }
    void AddGoodLumiRange(int run, int l0, int l1);

    // public configuration
    void SetFlag(TString flag, bool b=true) { flags[flag]=b; }
    bool isData=false;                                                 // to do gen matching, etc
    int firstEvent=-1;
    int lastEvent=-1;                                                    // max events to process; -1=>all
    ProcessType processType=kNone;                         // determine what to do the jet matching to

private:
    enum CorrectionType { //!< enum listing relevant corrections applied to MC
        cPU=0,         //!< true pu weight
        cPUUp,         //!< true pu weight Up
        cPUDown,       //!< true pu weight Down
	ZHEwkCorr,     //!< ZH Ewk Corr weight  
	ZHEwkCorrUp,   //!< ZH Ewk Corr weight Up  
	ZHEwkCorrDown, //!< ZH Ewk Corr weight Down  
	cLooseMuonId,
	cMediumMuonId,
	cTightMuonId,
	cLooseMuonIso, 
	cMediumMuonIso,
	cTightMuonIso,
	cTrackingMuon,
	cLooseElectronId, 
	cMediumElectronId,
	cTightElectronId,
	cTrackingElectron,
        cN
    };

    enum BTagType {
        bJetL=0,
        bJetM,
        bN
    };

    class btagcand {
        public:
            btagcand(unsigned int i, int f,double e,double cent,double up,double down) {
                idx = i;
                flav = f;
                eff = e;
                sf = cent;
                sfup = up;
                sfdown = down;
            }
            ~btagcand() { }
            int flav, idx;
            double eff, sf, sfup, sfdown;
    };

    bool PassGoodLumis(int run, int lumi);
    bool PassPreselection();
    void CalcBJetSFs(BTagType bt, int flavor, double eta, double pt, 
                         double eff, double uncFactor, double &sf, double &sfUp, double &sfDown);
    void EvalBTagSF(std::vector<btagcand> &cands, std::vector<double> &sfs,
                    GeneralLeptonicTree::BTagShift shift,GeneralLeptonicTree::BTagJet jettype, bool do2=false);
    void OpenCorrection(CorrectionType,TString,TString,int);
    double GetCorr(CorrectionType ct,double x, double y=0);
    void RegisterTrigger(TString path, std::vector<unsigned> &idxs); 

    int DEBUG = 0; //!< debug verbosity level
    std::map<TString,bool> flags;

    std::map<panda::GenParticle const*,float> genObjects;                 //!< particles we want to match the jets to, and the 'size' of the daughters
    panda::GenParticle const* MatchToGen(double eta, double phi, double r2, int pdgid=0);        //!< private function to match a jet; returns NULL if not found
    std::map<int,std::vector<LumiRange>> goodLumis;
    std::vector<panda::Particle*> matchPhos, matchEles, matchLeps;
    
    // CMSSW-provided utilities

    BTagCalibration *btagCalib=0;

    std::vector<BTagCalibrationReader*> btagReaders = std::vector<BTagCalibrationReader*>(bN,0); //!< maps BTagType to a reader 
    
    std::map<TString,JetCorrectionUncertainty*> ak8UncReader; //!< calculate JES unc on the fly
    JERReader *ak8JERReader; //!< fatjet jet energy resolution reader
    std::map<TString,JetCorrectionUncertainty*> ak4UncReader; //!< calculate JES unc on the fly
    std::map<TString,FactorizedJetCorrector*> ak4ScaleReader; //!< calculate JES on the fly
    JERReader *ak4JERReader; //!< fatjet jet energy resolution reader
    EraHandler eras = EraHandler(2016); //!< determining data-taking era, to be used for era-dependent JEC

    // files and histograms containing weights
    std::vector<TFile*> fCorrs = std::vector<TFile*>(cN,0); //!< files containing corrections
    std::vector<THCorr1*> h1Corrs = std::vector<THCorr1*>(cN,0); //!< histograms for binned corrections
    std::vector<THCorr2*> h2Corrs = std::vector<THCorr2*>(cN,0); //!< histograms for binned corrections

    // IO for the analyzer
    TFile *fOut;     // output file is owned by PandaLeptonicAnalyzer
    TTree *tOut;
    GeneralLeptonicTree *gt; // essentially a wrapper around tOut
    TH1F *hDTotalMCWeight=0;
    TH1D *hDDilPt;
    TTree *tIn=0;    // input tree to read
    unsigned int preselBits=0;

    // objects to read from the tree
    panda::Event event;

    // any extra signal weights we want
    std::vector<TString> wIDs;

};

#endif
