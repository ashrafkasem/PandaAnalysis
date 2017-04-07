#ifndef BTAGTREE_H
#define BTAGTREE_H

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "genericTree.h"
class BTagTree : public genericTree {
  public:
    BTagTree();
    ~BTagTree();
    void ReadTree(TTree *t);
    void WriteTree(TTree *t);
    void Reset();
    
	int runNumber=0;
	int lumiNumber=0;
	ULong64_t eventNumber=0;
	int idx=0;
	int npv=0;
	float weight=0;
	float pt=0;
	float eta=0;
	float genpt=0;
	int flavor=0;
	float csv=0;
//ENDDEF
};

#endif
