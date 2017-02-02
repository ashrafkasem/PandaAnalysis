#include "../interface/BTagTreeBuilder.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define DEBUG 0
using namespace std;

BTagTreeBuilder::BTagTreeBuilder() {
	bt = new BTagTree();
}

BTagTreeBuilder::~BTagTreeBuilder() {
	delete bt;
}

void BTagTreeBuilder::ResetBranches() {
	bt->Reset();
}

void BTagTreeBuilder::SetOutputPath(TString fOutName) {
	fOut = new TFile(fOutName,"RECREATE");
	tOut = new TTree("jets","jets");

	bt->Reset();
	bt->WriteTree(tOut);

}

void BTagTreeBuilder::SetInput(TString fInName) {
	fIn = new TFile(fInName);
	tIn = (TTree*)fIn->Get("events");
}

void BTagTreeBuilder::Terminate() {
	fOut->WriteTObject(tOut);
	fOut->Close();

	fIn->Close();
}

// run
void BTagTreeBuilder::Run() {
	
	ULong64_t eventNumber;
	int npv,	runNumber, lumiNumber;
	float mcWeight;
	float jet1Pt, jet1Eta, jet1GenPt, jet1CSV; 
	float jet1Flav; // will be int soon
	float jet2Pt, jet2Eta, jet2GenPt, jet2CSV; 
	float jet2Flav;

	tIn->SetBranchStatus("*",0);
	std::vector<TString> toRead = {"eventNumber","runNumber","lumiNumber","normalizedWeight","npv",
		                             "jet1Pt","jet1Eta","jet1GenPt","jet1CSV","jet1Flav",
		                             "jet2Pt","jet2Eta","jet2GenPt","jet2CSV","jet2Flav"};
	for (auto &b : toRead) 
		tIn->SetBranchStatus(b,1);

	tIn->SetBranchAddress("eventNumber",&eventNumber);
	tIn->SetBranchAddress("runNumber",&runNumber);
	tIn->SetBranchAddress("lumiNumber",&lumiNumber);
	tIn->SetBranchAddress("npv",&npv);
	tIn->SetBranchAddress("normalizedWeight",&mcWeight);
	tIn->SetBranchAddress("jet1Pt",&jet1Pt);
	tIn->SetBranchAddress("jet1Eta",&jet1Eta);
	tIn->SetBranchAddress("jet1GenPt",&jet1GenPt);
	tIn->SetBranchAddress("jet1CSV",&jet1CSV);
	tIn->SetBranchAddress("jet1Flav",&jet1Flav);
	tIn->SetBranchAddress("jet2Pt",&jet2Pt);
	tIn->SetBranchAddress("jet2Eta",&jet2Eta);
	tIn->SetBranchAddress("jet2GenPt",&jet2GenPt);
	tIn->SetBranchAddress("jet2CSV",&jet2CSV);
	tIn->SetBranchAddress("jet2Flav",&jet2Flav);

	unsigned nEvents = tIn->GetEntries();
	unsigned nZero = 0;
	if (lastEvent>=0 && lastEvent<(int)nEvents)
		nEvents = lastEvent;
	if (firstEvent>=0)
		nZero = firstEvent;

	unsigned int iE=nZero;
	ProgressReporter pr("BTagTreeBuilder::Run",&iE,&nEvents,10);
	TimeReporter tr("BTagTreeBuilder::Run",DEBUG);

	for (iE=nZero; iE!=nEvents; ++iE) {
		pr.Report();
		ResetBranches();

		tIn->GetEntry(iE);

		bt->eventNumber = eventNumber;
		bt->runNumber   = runNumber;
		bt->lumiNumber  = lumiNumber;
		bt->npv         = npv;
		bt->weight      = mcWeight;

		if (jet1Pt>20) {
			bt->pt         = jet1Pt;
			bt->eta        = jet1Eta;
			bt->genpt      = jet1GenPt;
			bt->flavor     = int(jet1Flav);
			bt->csv        = jet1CSV;
			bt->idx        = 1;
			tOut->Fill();
		}
		if (jet2Pt>20) {
			bt->pt         = jet2Pt;
			bt->eta        = jet2Eta;
			bt->genpt      = jet2GenPt;
			bt->flavor     = int(jet2Flav);
			bt->csv        = jet2CSV;
			bt->idx        = 2;
			tOut->Fill();
		}

	}
}

