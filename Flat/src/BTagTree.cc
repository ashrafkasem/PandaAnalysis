#include "../interface/BTagTree.h"

BTagTree::BTagTree() {
    	runNumber=0;
	lumiNumber=0;
	eventNumber=0;
	idx=0;
	npv=0;
	weight=0;
	pt=0;
	eta=0;
	genpt=0;
	flavor=0;
	csv=0;
//ENDCONST
}

BTagTree::~BTagTree() {
    //ENDDEST
}

void BTagTree::Reset() {
    	runNumber = 0;
	lumiNumber = 0;
	eventNumber = 0;
	idx = 0;
	npv = 0;
	weight = -1;
	pt = -1;
	eta = -1;
	genpt = -1;
	flavor = 0;
	csv = -1;
//ENDRESET
}

void BTagTree::ReadTree(TTree *t) {
      treePtr = t;
      treePtr->SetBranchStatus("*",0);
    	treePtr->SetBranchStatus("runNumber",1);
	treePtr->SetBranchAddress("runNumber",&runNumber);
	treePtr->SetBranchStatus("lumiNumber",1);
	treePtr->SetBranchAddress("lumiNumber",&lumiNumber);
	treePtr->SetBranchStatus("eventNumber",1);
	treePtr->SetBranchAddress("eventNumber",&eventNumber);
	treePtr->SetBranchStatus("idx",1);
	treePtr->SetBranchAddress("idx",&idx);
	treePtr->SetBranchStatus("npv",1);
	treePtr->SetBranchAddress("npv",&npv);
	treePtr->SetBranchStatus("weight",1);
	treePtr->SetBranchAddress("weight",&weight);
	treePtr->SetBranchStatus("pt",1);
	treePtr->SetBranchAddress("pt",&pt);
	treePtr->SetBranchStatus("eta",1);
	treePtr->SetBranchAddress("eta",&eta);
	treePtr->SetBranchStatus("genpt",1);
	treePtr->SetBranchAddress("genpt",&genpt);
	treePtr->SetBranchStatus("flavor",1);
	treePtr->SetBranchAddress("flavor",&flavor);
	treePtr->SetBranchStatus("csv",1);
	treePtr->SetBranchAddress("csv",&csv);
//ENDREAD
}

void BTagTree::WriteTree(TTree *t) {
      treePtr = t;
    	treePtr->Branch("runNumber",&runNumber,"runNumber/I");
	treePtr->Branch("lumiNumber",&lumiNumber,"lumiNumber/I");
	treePtr->Branch("eventNumber",&eventNumber,"eventNumber/l");
	treePtr->Branch("idx",&idx,"idx/I");
	treePtr->Branch("npv",&npv,"npv/I");
	treePtr->Branch("weight",&weight,"weight/F");
	treePtr->Branch("pt",&pt,"pt/F");
	treePtr->Branch("eta",&eta,"eta/F");
	treePtr->Branch("genpt",&genpt,"genpt/F");
	treePtr->Branch("flavor",&flavor,"flavor/I");
	treePtr->Branch("csv",&csv,"csv/F");
//ENDWRITE
}

