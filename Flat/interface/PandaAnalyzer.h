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
     kMonotop    =(1<<0),
     kMonohiggs  =(1<<2),
     kMonojet    =(1<<3),
     kTriggers   =(1<<4),
     kVBF        =(1<<5),
     kRecoil     =(1<<6)
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
        kSingleMuTrig     =(1<<3)
    };

    //////////////////////////////////////////////////////////////////////////////////////

    PandaAnalyzer();
    ~PandaAnalyzer();
    void Init(TTree *tree, TH1D *hweights);
    void SetOutputFile(TString fOutName);
    void ResetBranches();
    void Run();
    void Terminate();
    void SetDataDir(const char *s);
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

    std::map<TString,bool> flags;

    std::map<panda::GenParticle const*,float> genObjects;                 //!< particles we want to match the jets to, and the 'size' of the daughters
    panda::GenParticle const* MatchToGen(double eta, double phi, double r2, int pdgid=0);        //!< private function to match a jet; returns NULL if not found
    std::map<int,std::vector<LumiRange>> goodLumis;
    bool PassGoodLumis(int run, int lumi);
    bool PassPreselection();
    float getMSDcorr(Float_t puppipt, Float_t puppieta);
    std::vector<panda::Particle*> matchPhos, matchEles, matchLeps;
    
    // CMSSW-provided utilities

    void calcBJetSFs(TString readername, int flavor, double eta, double pt, 
                         double eff, double uncFactor, double &sf, double &sfUp, double &sfDown);
    BTagCalibration *btagCalib=0;
    BTagCalibration *sj_btagCalib=0;

    std::map<TString,BTagCalibrationReader*> btagReaders; //!< maps "JETTYPE_WP" to a reader 
                                                                                                                // I think we can load multiple flavors in a single reader 
    
    std::map<TString,JetCorrectionUncertainty*> ak8UncReader;            //!< calculate JES unc on the fly
    JERReader *ak8JERReader; //!< fatjet jet energy resolution reader
    EraHandler eras = EraHandler(2016); //!< determining data-taking era, to be used for era-dependent JEC

    std::vector<TFile*> openFiles; //!< anything that should be closed
    std::vector<void*> gc; //!< used for misc garbage collection

    // files and histograms containing weights
    TFile *fLepSF=0, *fLepRecoSF=0;
    TFile *fEleSF=0;
    THCorr2 *hEleVeto, *hEleTight;
    THCorr2 *hMuLooseLoPU, *hMuTightLoPU;
    THCorr2 *hMuLooseHiPU, *hMuTightHiPU;
    THCorr2 *hRecoEle;
    THCorr2 *hRecoMuLoPU,    *hRecoMuHiPU;
    TFile *fPhoSF=0;
    THCorr2 *hPho=0;

    TFile *fPU=0;
    THCorr1 *hPUWeight;

    TFile *fKFactor=0;
    THCorr1 *hZNLO, *hANLO, *hWNLO;
    THCorr1 *hZLO,    *hALO,    *hWLO;
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
    TFile *fOut;     // output file is owned by PandaAnalyzer
    TTree *tOut;
    GeneralTree *gt; // essentially a wrapper around tOut
    TTree *tIn=0;    // input tree to read
    unsigned int preselBits=0;

    // objects to read from the tree
    panda::Event event;

    // configuration read from output tree
    std::vector<double> betas;
    std::vector<int> Ns; 
    std::vector<int> orders;

};


#endif

