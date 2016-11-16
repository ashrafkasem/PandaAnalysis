#include "../interface/SFTreeBuilder.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define DEBUG 0
using namespace std;

SFTreeBuilder::SFTreeBuilder() {
}

SFTreeBuilder::~SFTreeBuilder() {
  delete weightTF;
  delete cutTF;
  for (auto *t : taggerTFs)
    delete t;
  delete taggerVals;
}

void SFTreeBuilder::ResetBranches() {
  runNumber=0;
  lumiNumber=0;
  eventNumber=0;
  weight=0;
  mSD=-1;
}

void SFTreeBuilder::SetOutputFile(TString fOutName, TString tOutName) {
  fOut = new TFile(fOutName,"RECREATE");
  tOut = new TTree(tOutName,tOutName);

  tOut->Branch("eventNumber",&eventNumber,"eventNumber/l");
  tOut->Branch("runNumber",&runNumber,"runNumber/I");
  tOut->Branch("lumiNumber",&lumiNumber,"lumiNumber/I");
  tOut->Branch("npv",&npv,"npv/I");
  tOut->Branch("pid",&pid,"pid/I");
  tOut->Branch("pt",&pt,"pt/F");
  tOut->Branch("eta",&eta,"eta/F");
  tOut->Branch("mSD",&mSD,"mSD/F");
  tOut->Branch("weight",&weight,"weight/F");
  tOut->Branch("pt",&pt,"pt/F");

  // by this point, all taggers should have been set
  taggerVals = new float[tagger_inNames.size()];
  for (unsigned iT=0; iT!=tagger_inNames.size(); ++iT) {
    TString outname = tagger_outNames[iT];
    tOut->Branch(outname,&(taggerVals[iT]),outname+"/F");
  }
}

void SFTreeBuilder::AddInputFile(TString fInName) {
  TFile *fIn = TFile::Open(fInName);
  TTree *t = (TTree*)fIn->FindObjectAny(inTreeName);
  processes.push_back(t);
}

void SFTreeBuilder::AddTagger(TString inName, TString outName) {
  if (outName=="")
    outName = inName;
  tagger_inNames.push_back(inName);
  tagger_outNames.push_back(outName);
}

void SFTreeBuilder::Terminate() {
  fOut->WriteTObject(tOut);
  fOut->Close();
}

// run
void SFTreeBuilder::Run() {

  for (TTree *tIn : processes) {
    pid++;

    weightTF = new TTreeFormula("weight",weightFormula,tIn); weightTF->SetQuickLoad(true);
    cutTF = new TTreeFormula("cut",cutFormula,tIn); cutTF->SetQuickLoad(true);
    unsigned nTaggers = tagger_inNames.size();
    for (unsigned iT=0; iT!=nTaggers; ++iT) {
      TTreeFormula *this_tf 
            = new TTreeFormula(tagger_outNames[iT], tagger_inNames[iT], tIn);
      this_tf->SetQuickLoad(true);
      taggerTFs.push_back(this_tf);
    }

    tIn->SetBranchAddress("eventNumber",&eventNumber);
    tIn->SetBranchAddress("runNumber",&runNumber);
    tIn->SetBranchAddress("lumiNumber",&lumiNumber);
    tIn->SetBranchAddress("npv",&npv);
    tIn->SetBranchAddress("fj1Pt",&pt);
    tIn->SetBranchAddress("fj1Eta",&eta);
    tIn->SetBranchAddress("fj1MSD",&mSD);

    // INITIALIZE --------------------------------------------------------------------------  
    unsigned int nEvents = tIn->GetEntries();
    unsigned int nZero = 0;
    if (lastEvent>=0 && lastEvent<(int)nEvents)
      nEvents = lastEvent;
    if (firstEvent>=0)
      nZero = firstEvent;

    // set up reporters
    unsigned int iE=0;
    ProgressReporter pr("SFTreeBuilder::Run",&iE,&nEvents,10);
    TimeReporter tr("SFTreeBuilder::Run",DEBUG);

    // EVENTLOOP --------------------------------------------------------------------------
    for (iE=nZero; iE!=nEvents; ++iE) {
      tr.Start();
      pr.Report();
      ResetBranches();
      tIn->GetEntry(iE);
      tr.TriggerEvent("GetEntry");

      // first check the cut
      if (!(cutTF->EvalInstance()))
        continue;

      tr.TriggerEvent("cut");

      weight = weightTF->EvalInstance();

      tr.TriggerEvent("weight");

      for (unsigned iT=0; iT!=nTaggers; ++iT) {
        double val = taggerTFs[iT]->EvalInstance();
        taggerVals[iT] = clean(val);
      }

      tr.TriggerEvent("taggers");

      tOut->Fill();
    }

    for (auto *t : taggerTFs)
      delete t;
    taggerTFs.clear();
    delete cutTF; cutTF=0;
    delete weightTF; weightTF=0;
  }

} // Run()

