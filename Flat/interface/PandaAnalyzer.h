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
     kRecoil     =(1<<6),
     kFatjet     =(1<<7)
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
        kSingleEleTrig =(1<<1),
        kSinglePhoTrig =(1<<2),
        kSingleMuTrig     =(1<<3)
    };

    //////////////////////////////////////////////////////////////////////////////////////

    PandaAnalyzer(int debug_=0);
    ~PandaAnalyzer();
    int Init(TTree *tree, TH1D *hweights, TTree *weightNames=0);
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
    enum CorrectionType { //!< enum listing relevant corrections applied to MC
        cNPV=0,       //!< npv weight
        cPU,          //!< true pu weight
        cEleVeto,     //!< monojet SF, Veto ID for e
        cEleTight,    //!< monojet SF, Tight ID for e
        cEleReco,     //!< monojet SF, tracking for e
        cMuLooseID,   //!< MUO POG SF, Loose ID for mu 
        cMuTightID,   //!< MUO POG SF, Tight ID for mu 
        cMuLooseIso,  //!< MUO POG SF, Loose Iso for mu 
        cMuTightIso,  //!< MUO POG SF, Tight Iso for mu 
        cMuReco,      //!< MUO POG SF, tracking for mu
        cPho,         //!< EGM POG SF, contains ID for gamma
        cTrigMET,     //!< MET trigger eff        
        cTrigMETZmm,  //!< Zmumu MET trigger eff
        cTrigEle,     //!< Ele trigger eff        
        cTrigPho,     //!< Pho trigger eff        
        cZNLO,        //!< NLO weights for QCD Z,W,A,A+2j
        cWNLO,
        cANLO,
        cANLO2j,
        cZEWK,        //!< EWK weights for QCD Z,W,A,A+2j
        cWEWK,
        cAEWK,
        cVBF_ZNLO,    //!< NLO weights for QCD Z,W in VBF phase space
        cVBF_WNLO,
        cVBFTight_ZNLO,    //!< NLO weights for QCD Z,W in tight VBF phase space
        cVBFTight_WNLO,
        cVBF_EWKZ,    //!< k-factors for EWK Z,W in VBF phase space
        cVBF_EWKW,
        cVBF_TrigMET, //!< MET trigger eff as a f'n of mjj/met 
        cVBF_TrigMETZmm,
        cBadECALJets,  //!< bad ECAL clusters to filter jets
        cN
    };

    enum BTagType {
        bJetL=0,
        bSubJetL,
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
    float GetMSDCorr(Float_t puppipt, Float_t puppieta);
    void CalcBJetSFs(BTagType bt, int flavor, double eta, double pt, 
                         double eff, double uncFactor, double &sf, double &sfUp, double &sfDown);
    void EvalBTagSF(std::vector<btagcand> &cands, std::vector<double> &sfs,
                    GeneralTree::BTagShift shift,GeneralTree::BTagJet jettype, bool do2=false);
    void OpenCorrection(CorrectionType,TString,TString,int);
    double GetCorr(CorrectionType ct,double x, double y=0);
    void RegisterTrigger(TString path, std::vector<unsigned> &idxs); 

    int DEBUG = 0; //!< debug verbosity level
    std::map<TString,bool> flags;

    std::map<panda::GenParticle const*,float> genObjects;                 //!< particles we want to match the jets to, and the 'size' of the daughters
    panda::GenParticle const* MatchToGen(double eta, double phi, double r2, int pdgid=0);        //!< private function to match a jet; returns NULL if not found
    std::map<int,std::vector<LumiRange>> goodLumis;
    std::vector<panda::Particle*> matchPhos, matchEles, matchLeps;
    
    // fastjet reclustering
    fastjet::JetDefinition *jetDef=0;
    fastjet::contrib::SoftDrop *softDrop=0;
    fastjet::AreaDefinition *areaDef=0;
    fastjet::GhostedAreaSpec *activeArea=0;

    // CMSSW-provided utilities

    BTagCalibration *btagCalib=0;
    BTagCalibration *sj_btagCalib=0;

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

    TFile *MSDcorr;
    TF1* puppisd_corrGEN;
    TF1* puppisd_corrRECO_cen;
    TF1* puppisd_corrRECO_for;

    // IO for the analyzer
    TFile *fOut;     // output file is owned by PandaAnalyzer
    TTree *tOut;
    GeneralTree *gt; // essentially a wrapper around tOut
    TH1F *hDTotalMCWeight=0;
    TTree *tIn=0;    // input tree to read
    unsigned int preselBits=0;

    // objects to read from the tree
    panda::Event event;

    // configuration read from output tree
    std::vector<int> ibetas;
    std::vector<int> Ns; 
    std::vector<int> orders;

    // any extra signal weights we want
    std::vector<TString> wIDs;

};


#endif

