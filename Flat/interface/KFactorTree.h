#ifndef KFACTORTREE_H
#define KFACTORTREE_H 

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "genericTree.h"
class KFactorTree : public genericTree {
  public:
    KFactorTree();
    ~KFactorTree();
    void ReadTree(TTree *t);
    void WriteTree(TTree *t);
    void Reset();
    
	ULong64_t eventNumber=0;
	float weight=0;
	int vid=0;
	float vpt=0;
	float veta=0;
	float vm=0;
	float ht=0;
	float lep1pt=0;
	float lep1eta=0;
	int lep1id=0;
	float lep2pt=0;
	float lep2eta=0;
	int lep2id=0;
	float jet1pt=0;
	float jet1eta=0;
	float jet2pt=0;
	float jet2eta=0;
	float mjj=0;
	float jjdeta=0;
	float jjdphi=0;
	int njet=0;
	float met=0;
//ENDDEF
};

#endif
