#ifndef GenAnalyzer_h
#define GenAnalyzer_h

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

#include "KFactorTree.h"

#include "PandaCore/Tools/interface/Common.h"


/////////////////////////////////////////////////////////////////////////////
// some misc definitions

/////////////////////////////////////////////////////////////////////////////
// GenAnalyzer definition
class GenAnalyzer {
public :
	// configuration enums
	enum Order {
	 kLO,
	 kNLO
	};

	enum Process { 
		kZ,
		kW
	};

	enum Plots {
		kVpt,
		kVeta,
		kMET,
		kMjj,
		kHT,
		kJet1pt,
		nPlots
	};

	//////////////////////////////////////////////////////////////////////////////////////

	GenAnalyzer();
	~GenAnalyzer();
	void AddInput(TTree *tree, double xsec, TString label);
	void SetOutputPath(TString fOutName, TString fOutFileName="skimmed.root"); //!< location of output file plus any generated histograms
	void ResetBranches();
	void Run();
	void Terminate();

	// public configuration
	void SetFlag(TString flag, bool b=true) { flags[flag]=b; }
	void SetCut(TString cut) { scut = cut; }
	int firstEvent=-1;
	int lastEvent=-1;													// max events to process; -1=>all
	Process processType;						 // determine which leps to look for
	Order order;												 // this determines the data format

private:
	struct InputTree {
		TTree *t;
		double xsec;
		TString label;
	};

	void RunLO(); 
	void RunNLO();

	std::map<TString,bool> flags;

	// IO for the analyzer
	TFile *fOut;	 // output file is owned by GenAnalyzer
	TTree *tOut;
	KFactorTree *kt; // essentially a wrapper around tOut
	std::vector<InputTree> tIns;	// input trees to read
	TString outpath;

	TString scut="1==1"; // turned into a TTreeFormula

	TH1D *hists[nPlots]{};

	unsigned int nEvents, nZero;
};


#endif

