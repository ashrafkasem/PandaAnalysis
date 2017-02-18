#include "../interface/KFactorTree.h"

KFactorTree::KFactorTree() {
    	eventNumber=0;
	weight=0;
	vid=0;
	vpt=0;
	veta=0;
	vm=0;
	ht=0;
	lep1pt=0;
	lep1eta=0;
	lep1id=0;
	lep2pt=0;
	lep2eta=0;
	lep2id=0;
	jet1pt=0;
	jet1eta=0;
	jet2pt=0;
	jet2eta=0;
	mjj=0;
	jjdeta=0;
	jjdphi=0;
	njet=0;
	met=0;
//ENDCONST
}

KFactorTree::~KFactorTree() {
    //ENDDEST
}

void KFactorTree::Reset() {
    	eventNumber = 0;
	weight = -1;
	vid = 0;
	vpt = -1;
	veta = -1;
	vm = -1;
	ht = -1;
	lep1pt = -1;
	lep1eta = -1;
	lep1id = 0;
	lep2pt = -1;
	lep2eta = -1;
	lep2id = 0;
	jet1pt = -1;
	jet1eta = -1;
	jet2pt = -1;
	jet2eta = -1;
	mjj = -1;
	jjdeta = -1;
	jjdphi = -1;
	njet = 0;
	met = -1;
//ENDRESET
}

void KFactorTree::ReadTree(TTree *t) {
      treePtr = t;
      treePtr->SetBranchStatus("*",0);
    	treePtr->SetBranchStatus("eventNumber",1);
	treePtr->SetBranchAddress("eventNumber",&eventNumber);
	treePtr->SetBranchStatus("weight",1);
	treePtr->SetBranchAddress("weight",&weight);
	treePtr->SetBranchStatus("vid",1);
	treePtr->SetBranchAddress("vid",&vid);
	treePtr->SetBranchStatus("vpt",1);
	treePtr->SetBranchAddress("vpt",&vpt);
	treePtr->SetBranchStatus("veta",1);
	treePtr->SetBranchAddress("veta",&veta);
	treePtr->SetBranchStatus("vm",1);
	treePtr->SetBranchAddress("vm",&vm);
	treePtr->SetBranchStatus("ht",1);
	treePtr->SetBranchAddress("ht",&ht);
	treePtr->SetBranchStatus("lep1pt",1);
	treePtr->SetBranchAddress("lep1pt",&lep1pt);
	treePtr->SetBranchStatus("lep1eta",1);
	treePtr->SetBranchAddress("lep1eta",&lep1eta);
	treePtr->SetBranchStatus("lep1id",1);
	treePtr->SetBranchAddress("lep1id",&lep1id);
	treePtr->SetBranchStatus("lep2pt",1);
	treePtr->SetBranchAddress("lep2pt",&lep2pt);
	treePtr->SetBranchStatus("lep2eta",1);
	treePtr->SetBranchAddress("lep2eta",&lep2eta);
	treePtr->SetBranchStatus("lep2id",1);
	treePtr->SetBranchAddress("lep2id",&lep2id);
	treePtr->SetBranchStatus("jet1pt",1);
	treePtr->SetBranchAddress("jet1pt",&jet1pt);
	treePtr->SetBranchStatus("jet1eta",1);
	treePtr->SetBranchAddress("jet1eta",&jet1eta);
	treePtr->SetBranchStatus("jet2pt",1);
	treePtr->SetBranchAddress("jet2pt",&jet2pt);
	treePtr->SetBranchStatus("jet2eta",1);
	treePtr->SetBranchAddress("jet2eta",&jet2eta);
	treePtr->SetBranchStatus("mjj",1);
	treePtr->SetBranchAddress("mjj",&mjj);
	treePtr->SetBranchStatus("jjdeta",1);
	treePtr->SetBranchAddress("jjdeta",&jjdeta);
	treePtr->SetBranchStatus("jjdphi",1);
	treePtr->SetBranchAddress("jjdphi",&jjdphi);
	treePtr->SetBranchStatus("njet",1);
	treePtr->SetBranchAddress("njet",&njet);
	treePtr->SetBranchStatus("met",1);
	treePtr->SetBranchAddress("met",&met);
//ENDREAD
}

void KFactorTree::WriteTree(TTree *t) {
      treePtr = t;
    	treePtr->Branch("eventNumber",&eventNumber,"eventNumber/l");
	treePtr->Branch("weight",&weight,"weight/F");
	treePtr->Branch("vid",&vid,"vid/I");
	treePtr->Branch("vpt",&vpt,"vpt/F");
	treePtr->Branch("veta",&veta,"veta/F");
	treePtr->Branch("vm",&vm,"vm/F");
	treePtr->Branch("ht",&ht,"ht/F");
	treePtr->Branch("lep1pt",&lep1pt,"lep1pt/F");
	treePtr->Branch("lep1eta",&lep1eta,"lep1eta/F");
	treePtr->Branch("lep1id",&lep1id,"lep1id/I");
	treePtr->Branch("lep2pt",&lep2pt,"lep2pt/F");
	treePtr->Branch("lep2eta",&lep2eta,"lep2eta/F");
	treePtr->Branch("lep2id",&lep2id,"lep2id/I");
	treePtr->Branch("jet1pt",&jet1pt,"jet1pt/F");
	treePtr->Branch("jet1eta",&jet1eta,"jet1eta/F");
	treePtr->Branch("jet2pt",&jet2pt,"jet2pt/F");
	treePtr->Branch("jet2eta",&jet2eta,"jet2eta/F");
	treePtr->Branch("mjj",&mjj,"mjj/F");
	treePtr->Branch("jjdeta",&jjdeta,"jjdeta/F");
	treePtr->Branch("jjdphi",&jjdphi,"jjdphi/F");
	treePtr->Branch("njet",&njet,"njet/I");
	treePtr->Branch("met",&met,"met/F");
//ENDWRITE
}

