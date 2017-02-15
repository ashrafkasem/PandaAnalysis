#ifndef BTagTreeBuilder_h
#define BTagTreeBuilder_h

// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>

#include "BTagTree.h"

#include "PandaCore/Tools/interface/Common.h"
#include "PandaCore/Tools/interface/TreeTools.h"


/////////////////////////////////////////////////////////////////////////////
// some misc definitions

/////////////////////////////////////////////////////////////////////////////
// BTagTreeBuilder definition
class BTagTreeBuilder {
public :

	BTagTreeBuilder();
	~BTagTreeBuilder();
	void SetInput(TString fInName);
	void SetOutputPath(TString fOutName); 
	void ResetBranches();
	void Run();
	void Terminate();

	// public configuration
	int firstEvent=-1;
	int lastEvent=-1;													// max events to process; -1=>all

private:
	// IO for the analyzer
	TFile *fOut=0;	 // output file is owned by BTagTreeBuilder
	TTree *tOut=0;
	BTagTree *bt=0; // essentially a wrapper around tOut
	TFile *fIn=0;
	TTree *tIn=0;	// input trees to read
};


#endif

