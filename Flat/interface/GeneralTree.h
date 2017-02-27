#ifndef GENERALTREE_H
#define GENERALTREE_H 

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TString.h"
#include "genericTree.h"
#include <map>

#define NJET 20
#define NSUBJET 2

class GeneralTree : public genericTree {
	private:
		std::vector<double> betas = {0.5, 1.0, 2.0, 4.0};
		std::vector<int> Ns = {1,2,3,4}; 
		std::vector<int> orders = {1,2,3};
		TString makeECFString(int order, int N, float beta) {
			return TString::Format("ECFN_%i_%i_%.2i",order,N,int(10*beta));
		}

	public:
		bool monohiggs=false, monojet=false, fatjet=true;
		GeneralTree();
		~GeneralTree();
		void ReadTree(TTree *t);
		void WriteTree(TTree *t);
		void Fill() { treePtr->Fill(); }
        void SetBranchStatus(const char *bname, bool status, UInt_t *ret=0) { treePtr->SetBranchStatus(bname,status,ret); }
		void Reset();

		std::vector<double> get_betas() const { return betas; }
		std::vector<int> get_Ns() const { return Ns; }
		std::vector<int> get_orders() const { return orders; }
		
	int runNumber=0;
	int lumiNumber=0;
	ULong64_t eventNumber=0;
	int npv=0;
	float mcWeight=0;
	float filter_maxRecoil=0;
	int filter_whichRecoil=-1;
	int trigger=0;
	int metFilter=0;
	float sf_ewkV=0;
	float sf_qcdV=0;
	float sf_qcdTT=0;
	float sf_lep=0;
	float sf_pho=0;
	float sf_lepReco=0;
	float sf_sjcsvWeightB=0;
	float sf_sjcsvWeightM=0;
	float sf_csvWeightB=0;
	float sf_csvWeightM=0;
	float sf_sjcsvWeightBUp=0;
	float sf_sjcsvWeightMUp=0;
	float sf_csvWeightBUp=0;
	float sf_csvWeightMUp=0;
	float sf_sjcsvWeightBDown=0;
	float sf_sjcsvWeightMDown=0;
	float sf_csvWeightBDown=0;
	float sf_csvWeightMDown=0;
	float sf_btag0=0;
	float sf_btag0BUp=0;
	float sf_btag0BDown=0;
	float sf_btag0MUp=0;
	float sf_btag0MDown=0;
	float sf_btag1=0;
	float sf_btag1BUp=0;
	float sf_btag1BDown=0;
	float sf_btag1MUp=0;
	float sf_btag1MDown=0;
	float sf_btagGT0=0;
	float sf_btagGT0BUp=0;
	float sf_btagGT0BDown=0;
	float sf_btagGT0MUp=0;
	float sf_btagGT0MDown=0;
	float sf_btag2=0;
	float sf_btag2BUp=0;
	float sf_btag2BDown=0;
	float sf_btag2MUp=0;
	float sf_btag2MDown=0;
	float sf_btag0_alt=0;
	float sf_btag0BUp_alt=0;
	float sf_btag0BDown_alt=0;
	float sf_btag0MUp_alt=0;
	float sf_btag0MDown_alt=0;
	float sf_btag1_alt=0;
	float sf_btag1BUp_alt=0;
	float sf_btag1BDown_alt=0;
	float sf_btag1MUp_alt=0;
	float sf_btag1MDown_alt=0;
	float sf_btagGT0_alt=0;
	float sf_btagGT0BUp_alt=0;
	float sf_btagGT0BDown_alt=0;
	float sf_btagGT0MUp_alt=0;
	float sf_btagGT0MDown_alt=0;
	float sf_btag2_alt=0;
	float sf_btag2BUp_alt=0;
	float sf_btag2BDown_alt=0;
	float sf_btag2MUp_alt=0;
	float sf_btag2MDown_alt=0;
	float sf_sjbtag0=0;
	float sf_sjbtag0BUp=0;
	float sf_sjbtag0BDown=0;
	float sf_sjbtag0MUp=0;
	float sf_sjbtag0MDown=0;
	float sf_sjbtag1=0;
	float sf_sjbtag1BUp=0;
	float sf_sjbtag1BDown=0;
	float sf_sjbtag1MUp=0;
	float sf_sjbtag1MDown=0;
	float sf_sjbtagGT0=0;
	float sf_sjbtagGT0BUp=0;
	float sf_sjbtagGT0BDown=0;
	float sf_sjbtagGT0MUp=0;
	float sf_sjbtagGT0MDown=0;
	float sf_sjbtag2=0;
	float sf_sjbtag2BUp=0;
	float sf_sjbtag2BDown=0;
	float sf_sjbtag2MUp=0;
	float sf_sjbtag2MDown=0;
	float sf_eleTrig=0;
	float sf_phoTrig=0;
	float sf_metTrig=0;
	float sf_pu=0;
	float sf_tt=0;
	float sf_tt_ext=0;
	float sf_tt_bound=0;
	float sf_tt8TeV=0;
	float sf_tt8TeV_ext=0;
	float sf_tt8TeV_bound=0;
	float sf_phoPurity=0;
	float finalWeight=0;
	float pfmet=0;
	float pfmetphi=0;
	float pfmetnomu=0;
	float puppimet=0;
	float puppimetphi=0;
	float calomet=0;
	float calometphi=0;
	float pfcalobalance=0;
	float sumET=0;
	float trkmet=0;
	float UWmag=0;
	float UWphi=0;
	float UZmag=0;
	float UZphi=0;
	float UAmag=0;
	float UAphi=0;
	float Uperp=0;
	float Upara=0;
	float pfUWmag=0;
	float pfUWphi=0;
	float pfUZmag=0;
	float pfUZphi=0;
	float pfUAmag=0;
	float pfUAphi=0;
	float pfUperp=0;
	float pfUpara=0;
	float dphipfmet=0;
	float dphipuppimet=0;
	float dphiUW=0;
	float dphiUZ=0;
	float dphiUA=0;
	float dphipfUW=0;
	float dphipfUZ=0;
	float dphipfUA=0;
	float trueGenBosonPt=0;
	float genBosonPt=0;
	float genBosonEta=0;
	float genBosonMass=0;
	float genBosonPhi=0;
	float genWPlusPt=0;
	float genWMinusPt=0;
	float genWPlusEta=0;
	float genWMinusEta=0;
	float genTopIsHad=0;
	float genTopPt=0;
	float genAntiTopIsHad=0;
	float genAntiTopPt=0;
	float genTopEta=0;
	float genAntiTopEta=0;
	float genTTPt=0;
	float genTTEta=0;
	int nJet=0;
	int nIsoJet=0;
	float jet1Phi=0;
	int jet1Flav=0;
	float jet1GenPt=0;
	float jet1Pt=0;
	float jet1Eta=0;
	float jet1CSV=0;
	float jet2Phi=0;
	int jet2Flav=0;
	float jet2GenPt=0;
	float jet2Pt=0;
	float jet2Eta=0;
	float jet2CSV=0;
	int jet1IsTight=0;
	float isojet1Pt=0;
	float isojet1CSV=0;
	int isojet1Flav=0;
	float isojet2Pt=0;
	float isojet2CSV=0;
	int isojet2Flav=0;
	float jet12Mass=0;
	float jet12DEta=0;
	int jetNBtags=0;
	int isojetNBtags=0;
	int nFatjet=0;
	float fj1Tau32=0;
	float fj1Tau21=0;
	float fj1Tau32SD=0;
	float fj1Tau21SD=0;
	float fj1MSD=0;
	float fj1MSD_corr=0;
	float fj1Pt=0;
	float fj1PtScaleUp=0;
	float fj1PtScaleDown=0;
	float fj1PtSmeared=0;
	float fj1PtSmearedUp=0;
	float fj1PtSmearedDown=0;
	float fj1MSDScaleUp=0;
	float fj1MSDScaleDown=0;
	float fj1MSDSmeared=0;
	float fj1MSDSmearedUp=0;
	float fj1MSDSmearedDown=0;
	float fj1Phi=0;
	float fj1Eta=0;
	float fj1M=0;
	float fj1MaxCSV=0;
	float fj1MinCSV=0;
	float fj1DoubleCSV=0;
	float fj1GenPt=0;
	float fj1GenSize=0;
	int fj1IsMatched=0;
	float fj1sjPt[NSUBJET];
	float fj1sjEta[NSUBJET];
	float fj1sjPhi[NSUBJET];
	float fj1sjM[NSUBJET];
	float fj1sjCSV[NSUBJET];
	float fj1sjQGL[NSUBJET];
	float fj1GenWPt=0;
	float fj1GenWSize=0;
	int fj1IsWMatched=0;
	int fj1HighestPtGen=0;
	float fj1HighestPtGenPt=-1;
	int fj1IsTight=0;
	int fj1IsLoose=0;
	float fj1RawPt=0;
	int fj1IsHF=0;
	float fj1HTTMass=0;
	float fj1HTTFRec=0;
	int fj1IsClean=0;
	int isHF=0;
	int nLoosePhoton=0;
	int nTightPhoton=0;
	int loosePho1IsTight=0;
	float loosePho1Pt=0;
	float loosePho1Eta=0;
	float loosePho1Phi=0;
	int nLooseLep=0;
	int nLooseElectron=0;
	int nLooseMuon=0;
	int nTightLep=0;
	int nTightElectron=0;
	int nTightMuon=0;
	int looseLep1PdgId=0;
	int looseLep2PdgId=0;
	int looseLep1IsTight=0;
	int looseLep2IsTight=0;
	int looseLep1IsHLTSafe=0;
	int looseLep2IsHLTSafe=0;
	float looseLep1Pt=0;
	float looseLep1Eta=0;
	float looseLep1Phi=0;
	float looseLep2Pt=0;
	float looseLep2Eta=0;
	float looseLep2Phi=0;
	float diLepMass=0;
	float mT=0;
	int nTau=0;
	float jetPt[NJET];
	float jetEta[NJET];
	float jetPhi[NJET];
	float jetE[NJET];
	float jetCSV[NJET];
	float jetIso[NJET];
	float jetQGL[NJET];
	float hbbpt;
	float hbbeta;
	float hbbphi;
	float hbbm;
	int hbbjtidx[2];
	float scaleUp;
	float scaleDown;
	float pdfUp;
	float pdfDown;
	float scale[6];
		
	std::map<TString,float> fj1ECFNs;

//ENDDEF
};

#endif
