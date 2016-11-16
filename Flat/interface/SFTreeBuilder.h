#ifndef SFTreeBuilder_h
#define SFTreeBuilder_h

// STL
#include "vector"
#include "map"
#include <string>
#include <cmath>

// ROOT
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TTreeFormula.h>

#include "PandaCore/Tools/interface/Common.h"
/////////////////////////////////////////////////////////////////////////////
// some misc definitions

/////////////////////////////////////////////////////////////////////////////
// SFTreeBuilder definition
class SFTreeBuilder {
public :
  //////////////////////////////////////////////////////////////////////////////////////

  SFTreeBuilder();
  ~SFTreeBuilder();
  void AddInputFile(TString fInName);
  void SetOutputFile(TString fOutName, TString tOutName);
  void AddTagger(TString inName, TString outName="");
  void Run();
  void Terminate();

  // public configuration
  TString cutFormula;
  TString weightFormula;
  TString inTreeName="events";
  int firstEvent=-1;
  int lastEvent=-1;                          // max events to process; -1=>all

private:

  void ResetBranches();

  // vector of input trees
  std::vector<TTree*> processes; 

  // taggers
  std::vector<TString> tagger_inNames, tagger_outNames;
  std::vector<TTreeFormula*> taggerTFs;
  float *taggerVals=0;

  TTreeFormula *weightTF=0;
  TTreeFormula *cutTF=0;

  ULong64_t eventNumber;
  int runNumber=0, lumiNumber=0, npv=0, pid=-1;
  float weight=0, pt=0, eta=0, mSD=0; 


  // IO for the analyzer
  TFile *fOut;   // output file is owned by SFTreeBuilder
  TTree *tOut;

};


#endif

