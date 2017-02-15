#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define DEBUG 0
using namespace panda;
using namespace std;

PandaAnalyzer::PandaAnalyzer() {
	gt = new GeneralTree();
	betas = gt->get_betas(); 
	Ns = gt->get_Ns(); 
	orders = gt->get_orders(); 
	flags["fatjet"]      = true;
	flags["puppi"]       = true;
	flags["monohiggs"]   = false;
	flags["monojet"]     = false;
	flags["firstGen"]    = true;
	flags["applyJSON"]   = true;
}


PandaAnalyzer::~PandaAnalyzer() {
}


void PandaAnalyzer::ResetBranches() {
	genObjects.clear();
	matchPhos.clear();
	matchEles.clear();
	matchLeps.clear();
	gt->Reset();
}


void PandaAnalyzer::SetOutputFile(TString fOutName) {
	fOut = new TFile(fOutName,"RECREATE");
	tOut = new TTree("events","events");

	gt->monohiggs = flags["monohiggs"];
	gt->monojet   = flags["monojet"];
	gt->fatjet    = flags["fatjet"];
	gt->Reset(); // to be extra safe and fill the map before setting addresses
	gt->WriteTree(tOut);

}


void PandaAnalyzer::Init(TTree *t, TTree *infotree)
{
	if (!t) return;
	tIn = t;
	t->SetBranchAddress("event",&event);
	TString jetname = (flags["puppi"]) ? "puppi" : "chs";
	if (flags["fatjet"])
		t->SetBranchAddress(jetname+"CA15",&fatjets);
	t->SetBranchAddress("chsAK4",&jets);
	//t->SetBranchAddress(jetname+"AK4",&jets);
	t->SetBranchAddress("electron",&electrons);
	t->SetBranchAddress("muon",&muons);
	t->SetBranchAddress("tau",&taus);
	t->SetBranchAddress("photon",&photons);
	t->SetBranchAddress("pfmet",&pfmet);
	t->SetBranchAddress("puppimet",&puppimet);
	if (!isData) {
		t->SetBranchAddress("gen",&genparts);	
		t->SetBranchAddress("geninfo",&geninfo);
	}
	TH1F *hDTotalMCWeight = new TH1F("hDTotalMCWeight","hDTotalMCWeight",4,-2,2);
	TString val("fabs(info.mcWeight)/(info.mcWeight)");
	if (isData)	val = "1";
	infotree->Draw(val+">>hDTotalMCWeight",val);
	fOut->WriteTObject(hDTotalMCWeight,"hDTotalMCWeight");
}


PGenParticle *PandaAnalyzer::MatchToGen(double eta, double phi, double radius, int pdgid) {
	PGenParticle *found=NULL;
	double r2 = radius*radius;
	pdgid = abs(pdgid);

	unsigned int counter=0;
	for (map<PGenParticle*,float>::iterator iG=genObjects.begin();
				iG!=genObjects.end(); ++iG) {
		if (found!=NULL)
			break;
		if (pdgid!=0 && abs(iG->first->pdgid)!=pdgid)
			continue;
		if (DeltaR2(eta,phi,iG->first->eta,iG->first->phi)<r2) 
			found = iG->first;
	}

	return found;
}


void PandaAnalyzer::Terminate() {
	fOut->WriteTObject(tOut);
	fOut->Close();

	for (TFile *f :openFiles)
		f->Close();
	openFiles.clear();

	delete btagCalib;
	delete btagReaders["jet_L"];

	delete btagReaders["jet_M"];

	delete sj_btagCalib;
	delete btagReaders["sj_L"];

	for (auto& iter : ak8UncReader) 
		delete iter.second;

	delete ak8JERReader;
}


void PandaAnalyzer::SetDataDir(const char *s) {
	TString dirPath(s);

	fEleTrigB		= new TFile(dirPath+"/trigger_eff/ele_trig_lowpt_rebinned.root");
	fEleTrigE		= new TFile(dirPath+"/trigger_eff/ele_trig_lowpt_rebinned.root");
	fPhoTrig		 = new TFile(dirPath+"/trigger_eff/pho_trig.root");
	fMetTrig		 = new TFile(dirPath+"/trigger_eff/met_trig.root");
	fEleTrigLow	= new TFile(dirPath+"/trigger_eff/ele_trig_lowpt.root");
	openFiles.push_back(fEleTrigB);
	openFiles.push_back(fEleTrigE);
	openFiles.push_back(fPhoTrig);
	openFiles.push_back(fMetTrig);
	openFiles.push_back(fEleTrigLow);

	fLepSF				= new TFile(dirPath+"/scalefactors_lepton_moriond.root");
	fLepRecoSF		= new TFile(dirPath+"/scalefactors_reco_lepton_moriond.root");
	fEleSF				= new TFile(dirPath+"/scalefactors_80x_egpog_37ifb.root");
	openFiles.push_back(fLepSF);
	openFiles.push_back(fLepRecoSF);
	openFiles.push_back(fEleSF);

	fPhoSF	 = new TFile(dirPath+"/scalefactors_80x_medium_photon_37ifb.root");
	openFiles.push_back(fPhoSF);

	fPU			= new TFile(dirPath+"/puWeights_npv.root");
	openFiles.push_back(fPU);

	fKFactor = new TFile(dirPath+"/kfactors.root");
	openFiles.push_back(fKFactor);

	hEleTrigB = new THCorr1((TH1D*) fEleTrigB->Get("h_num"));
	hEleTrigE = new THCorr1((TH1D*) fEleTrigE->Get("h_num_endcap"));
	gc.push_back(hEleTrigB); gc.push_back(hEleTrigE);

	hPhoTrig		= new THCorr1((TH1D*) fPhoTrig->Get("h_num"));
	hMetTrig		= new THCorr1((TH1D*) fMetTrig->Get("numer"));
	hEleTrigLow = new THCorr2((TH2D*) fEleTrigLow->Get("hEffEtaPt"));
	gc.push_back(hPhoTrig); gc.push_back(hMetTrig); gc.push_back(hEleTrigLow);

	hEleVeto    	= new THCorr2((TH2D*) fEleSF->Get("scalefactors_Veto_Electron"));
	hEleTight    	= new THCorr2((TH2D*) fEleSF->Get("scalefactors_Tight_Electron"));
	hMuLooseLoPU	= new THCorr2((TH2D*) fLepSF->Get("scaleFactor_muon_looseid_pu_0_17"));
	hMuLooseHiPU	= new THCorr2((TH2D*) fLepSF->Get("scaleFactor_muon_looseid_pu_17_50"));
	hMuTightLoPU	= new THCorr2((TH2D*) fLepSF->Get("scaleFactor_muon_tightid_pu_0_17"));
	hMuTightHiPU	= new THCorr2((TH2D*) fLepSF->Get("scaleFactor_muon_tightid_pu_17_50"));

	hRecoEle  	= new THCorr2((TH2D*) fEleSF->Get("scalefactors_Reco_Electron"));
	hRecoMuLoPU	= new THCorr2((TH2D*) fLepRecoSF->Get("scaleFactor_muon_trackerid_pu_0_17"));
	hRecoMuHiPU	= new THCorr2((TH2D*) fLepRecoSF->Get("scaleFactor_muon_trackerid_pu_17_50"));

	hPho = new THCorr2((TH2D*) fPhoSF->Get("EGamma_SF2D"));
	gc.push_back(hPho);

	hPUWeight = new THCorr1((TH1D*)fPU->Get("data_npv_Wmn"));
	gc.push_back(hPUWeight);

	hZNLO = new THCorr1((TH1D*)fKFactor->Get("ZJets_012j_NLO/nominal"));
	hWNLO = new THCorr1((TH1D*)fKFactor->Get("WJets_012j_NLO/nominal"));
	hANLO = new THCorr1((TH1D*)fKFactor->Get("GJets_1j_NLO/nominal_G"));

	hZLO	= new THCorr1((TH1D*)fKFactor->Get("ZJets_LO/inv_pt"));
	hWLO	= new THCorr1((TH1D*)fKFactor->Get("WJets_LO/inv_pt"));
	hALO	= new THCorr1((TH1D*)fKFactor->Get("GJets_LO/inv_pt_G"));
 
	hZEWK = new THCorr1((TH1D*)fKFactor->Get("EWKcorr/Z"));
	hWEWK = new THCorr1((TH1D*)fKFactor->Get("EWKcorr/W"));
	hAEWK = new THCorr1((TH1D*)fKFactor->Get("EWKcorr/photon"));

	hZEWK->GetHist()->Divide(hZNLO->GetHist());	 
	hWEWK->GetHist()->Divide(hWNLO->GetHist());	 
	hAEWK->GetHist()->Divide(hANLO->GetHist());
	hZNLO->GetHist()->Divide(hZLO->GetHist());		
	hWNLO->GetHist()->Divide(hWLO->GetHist());		
	hANLO->GetHist()->Divide(hALO->GetHist());

	fCSVLF = new TFile(dirPath+"/csvWeights/csvweight_fake.root"); openFiles.push_back(fCSVLF);
	hCSVLF = new THCorr1( (TH1D*)fCSVLF->Get("hratio") ); gc.push_back(hCSVLF);

	fCSVHF = new TFile(dirPath+"/csvWeights/csvweight_tag_iterative.root"); openFiles.push_back(fCSVHF);
	hCSVHF = new THCorr1( (TH1D*)fCSVHF->Get("hratio") ); gc.push_back(hCSVHF);

	btagCalib = new BTagCalibration("csvv2",(dirPath+"/CSVv2_Moriond17_B_H.csv").Data());
	btagReaders["jet_L"] = new BTagCalibrationReader(BTagEntry::OP_LOOSE,"central",{"up","down"});
	btagReaders["jet_L"]->load(*btagCalib,BTagEntry::FLAV_B,"comb");
	btagReaders["jet_L"]->load(*btagCalib,BTagEntry::FLAV_C,"comb");
	btagReaders["jet_L"]->load(*btagCalib,BTagEntry::FLAV_UDSG,"incl");

	sj_btagCalib = new BTagCalibration("csvv2",(dirPath+"/subjet_CSVv2_Moriond17_B_H.csv").Data());
	btagReaders["sj_L"] = new BTagCalibrationReader(BTagEntry::OP_LOOSE,"central",{"up","down"});
	btagReaders["sj_L"]->load(*sj_btagCalib,BTagEntry::FLAV_B,"lt");
	btagReaders["sj_L"]->load(*sj_btagCalib,BTagEntry::FLAV_C,"lt");
	btagReaders["sj_L"]->load(*sj_btagCalib,BTagEntry::FLAV_UDSG,"incl");

	btagReaders["jet_M"] = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"});
	btagReaders["jet_M"]->load(*btagCalib,BTagEntry::FLAV_B,"comb");
	btagReaders["jet_M"]->load(*btagCalib,BTagEntry::FLAV_C,"comb");
	btagReaders["jet_M"]->load(*btagCalib,BTagEntry::FLAV_UDSG,"incl");

	if (flags["monohiggs"]) {
		MSDcorr = new TFile(dirPath+"/puppiCorr.root");
		puppisd_corrGEN = (TF1*)MSDcorr->Get("puppiJECcorr_gen");;
		puppisd_corrRECO_cen = (TF1*)MSDcorr->Get("puppiJECcorr_reco_0eta1v3");
		puppisd_corrRECO_for = (TF1*)MSDcorr->Get("puppiJECcorr_reco_1v3eta2v5");
	}

	ak8UncReader["MC"] = new JetCorrectionUncertainty(
			  (dirPath+"/jec/23Sep2016V2/Spring16_23Sep2016V2_MC_Uncertainty_AK8PFPuppi.txt").Data()
			);
	std::vector<TString> eraGroups = {"BCD","EF","G","H"};
	for (auto e : eraGroups) {
		ak8UncReader["data"+e] = new JetCorrectionUncertainty(
				  (dirPath+"/jec/23Sep2016V2/Spring16_23Sep2016"+e+"V2_DATA_Uncertainty_AK8PFPuppi.txt").Data()
				);
	}

	ak8JERReader = new JERReader(dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_SF_AK8PFPuppi.txt",
		                           dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt");

	// load only L2L3 JEC
	/*
	std::string jecPath = (dirPath+"/jec/").Data();
	std::vector<JetCorrectorParameters> mcParams;
	std::vector<JetCorrectorParameters> dataParams;
	mcParams.push_back(JetCorrectorParameters(jecPath + "Spring16_25nsV6_MC_L2Relative_AK8PFPuppi.txt"));
	mcParams.push_back(JetCorrectorParameters(jecPath + "Spring16_25nsV6_MC_L3Absolute_AK8PFPuppi.txt"));
	mcParams.push_back(JetCorrectorParameters(jecPath + "Spring16_25nsV6_MC_L2L3Residual_AK8PFPuppi.txt"));
	dataParams.push_back(JetCorrectorParameters(jecPath + "Spring16_25nsV6_DATA_L2Relative_AK8PFPuppi.txt"));
	dataParams.push_back(JetCorrectorParameters(jecPath + "Spring16_25nsV6_DATA_L3Absolute_AK8PFPuppi.txt"));
	dataParams.push_back(JetCorrectorParameters(jecPath + "Spring16_25nsV6_DATA_L2L3Residual_AK8PFPuppi.txt"));
	ak8MCCorrector = new FactorizedJetCorrector(mcParams);
	ak8DataCorrector = new FactorizedJetCorrector(dataParams);
	*/
//	ak8jec = new JetCorrectorParameters((dirPath+"/Spring16_25nsV6_MC_Uncertainty_AK8PFPuppi.txt").Data());
//	ak8UncReader = new JetCorrectionUncertainty(*ak8jec);
}


void PandaAnalyzer::AddGoodLumiRange(int run, int l0, int l1) {
	auto run_ = goodLumis.find(run); 
	if (run_==goodLumis.end()) { // don't know about this run yet
		std::vector<LumiRange> newLumiList;
		newLumiList.emplace_back(l0,l1);
		goodLumis[run] = newLumiList;
	} else {
		run_->second.emplace_back(l0,l1);
	}
}


bool PandaAnalyzer::PassGoodLumis(int run, int lumi) {
	auto run_ = goodLumis.find(run);
	if (run_==goodLumis.end()) {
		// matched no run
		if (DEBUG) PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i",run));
		return false;
	}

	// found the run, now look for a lumi range
	for (auto &range : run_->second) {
		if (range.Contains(lumi)) {
			if (DEBUG) PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Accepting run=%i, lumi=%i",run,lumi));
			return true;
		}
	}

	// matched no lumi range
	if (DEBUG) PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i, lumi=%i",run,lumi));
	return false;
}


bool PandaAnalyzer::PassPreselection() {
	if (preselBits==0)
		return true;
	bool isGood=false;

	if (preselBits & kRecoil) {
		if ( (gt->puppimet>200 || gt->UZmag>200 || gt->UWmag>200 || gt->UAmag>200) ||
					(gt->pfmet>200 || gt->pfUZmag>200 || gt->pfUWmag>200 || gt->pfUAmag>200) ) {
					isGood = true;
		}
	}
	if (preselBits & kMonotop) {
		if (gt->nFatjet>=1 && gt->fj1Pt>200) { 
			if ( (gt->puppimet>200 || gt->UZmag>200 || gt->UWmag>200 || gt->UAmag>200) ||
						(gt->pfmet>200 || gt->pfUZmag>200 || gt->pfUWmag>200 || gt->pfUAmag>200) ) {
						isGood = true;
			}
		}
	}
	if (preselBits & kMonojet) {
		if (true) { 
			if ( (gt->puppimet>200 || gt->UZmag>200 || gt->UWmag>200 || gt->UAmag>200) ||
						(gt->pfmet>200 || gt->pfUZmag>200 || gt->pfUWmag>200 || gt->pfUAmag>200) ) {
						isGood = true;
			}
		}
	}
	if (preselBits & kMonohiggs) {
		if ((gt->nFatjet>=1 && gt->fj1Pt>200) || gt->hbbpt>150 ) { 
			if ( (gt->puppimet>175 || gt->UZmag>175 || gt->UWmag>175 || gt->UAmag>175) ||
						(gt->pfmet>175 || gt->pfUZmag>175 || gt->pfUWmag>175 || gt->pfUAmag>175) ) {
						isGood = true;
			}
		}
	}
	
	return isGood;
}


void PandaAnalyzer::calcBJetSFs(TString readername, int flavor, 
																double eta, double pt, double eff, double uncFactor,
																double &sf, double &sfUp, double &sfDown) {
	if (flavor==5) {
		sf     = btagReaders[readername]->eval_auto_bounds("central",BTagEntry::FLAV_B,eta,pt);
		sfUp   = btagReaders[readername]->eval_auto_bounds("up",BTagEntry::FLAV_B,eta,pt);
		sfDown = btagReaders[readername]->eval_auto_bounds("down",BTagEntry::FLAV_B,eta,pt);
	} else if (flavor==4) {
		sf     = btagReaders[readername]->eval_auto_bounds("central",BTagEntry::FLAV_C,eta,pt);
		sfUp   = btagReaders[readername]->eval_auto_bounds("up",BTagEntry::FLAV_C,eta,pt);
		sfDown = btagReaders[readername]->eval_auto_bounds("down",BTagEntry::FLAV_C,eta,pt);
	} else {
		sf     = btagReaders[readername]->eval_auto_bounds("central",BTagEntry::FLAV_UDSG,eta,pt);
		sfUp   = btagReaders[readername]->eval_auto_bounds("up",BTagEntry::FLAV_UDSG,eta,pt);
		sfDown = btagReaders[readername]->eval_auto_bounds("down",BTagEntry::FLAV_UDSG,eta,pt);
	}

	sfUp = uncFactor*(sfUp-sf)+sf;
	sfDown = uncFactor*(sfDown-sf)+sf;
	return;
}


float PandaAnalyzer::getMSDcorr(Float_t puppipt, Float_t puppieta) { 

	float genCorr	= 1.;
	float recoCorr = 1.;
	float totalWeight = 1.;

	genCorr =	puppisd_corrGEN->Eval( puppipt );
	if( fabs(puppieta)	<= 1.3 ){
		recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
	}
	else{
		recoCorr = puppisd_corrRECO_for->Eval( puppipt );
	}
	totalWeight = genCorr * recoCorr;

	return totalWeight;
} 


// run
void PandaAnalyzer::Run() {

	// INITIALIZE --------------------------------------------------------------------------
	
	unsigned int nEvents = tIn->GetEntries();
	unsigned int nZero = 0;
	if (lastEvent>=0 && lastEvent<(int)nEvents)
		nEvents = lastEvent;
	if (firstEvent>=0)
		nZero = firstEvent;

	if (!fOut || !tIn) {
		PError("PandaAnalyzer::Run","NOT SETUP CORRECTLY");
		exit(1);
	}

	// get bounds
	float genBosonPtMin=150, genBosonPtMax=1000;
	if (!isData) {
		genBosonPtMin = hZNLO->GetHist()->GetBinCenter(1); 
		genBosonPtMax = hZNLO->GetHist()->GetBinCenter(hZNLO->GetHist()->GetNbinsX());
	}

	// these are bins of b-tagging eff in pT and eta, derived in 8024 TT MC
	std::vector<double> vbtagpt {20.0,50.0,80.0,120.0,200.0,300.0,400.0,500.0,700.0,1000.0};
	std::vector<double> vbtageta {0.0,0.5,1.5,2.5};
	std::vector<std::vector<double>> lfeff  = {{0.081,0.065,0.060,0.063,0.072,0.085,0.104,0.127,0.162},
                                         		 {0.116,0.097,0.092,0.099,0.112,0.138,0.166,0.185,0.222},
                                         		 {0.173,0.145,0.149,0.175,0.195,0.225,0.229,0.233,0.250}};
	std::vector<std::vector<double>> ceff = {{0.377,0.389,0.391,0.390,0.391,0.375,0.372,0.392,0.435},
                           		             {0.398,0.407,0.416,0.424,0.424,0.428,0.448,0.466,0.500},
                           		             {0.375,0.389,0.400,0.425,0.437,0.459,0.481,0.534,0.488}};
	std::vector<std::vector<double>> beff = {{0.791,0.815,0.825,0.835,0.821,0.799,0.784,0.767,0.760},
                          		             {0.794,0.816,0.829,0.836,0.823,0.804,0.798,0.792,0.789},
                          		             {0.739,0.767,0.780,0.789,0.776,0.771,0.779,0.787,0.806}};
	Binner btagpt(vbtagpt);
	Binner btageta(vbtageta);

	std::vector<unsigned> scaleIndices {1,2,3,4,6,8};

	// these are triggers. at some point these ought to be read from the file
	std::vector<unsigned int> metTriggers {0,1,2,3,4,5,6,7,8,9,10,11};
	std::vector<unsigned int> eleTriggers {27,28,29,30,31,32,33,34,35,36};
	std::vector<unsigned int> phoTriggers {37,38,39,40,41,42,43,44};
	std::vector<unsigned int> muTriggers  {13,14,15,16,17,18,19,20,21,22,23,24,25,26};

	// set up reporters
	unsigned int iE=0;
	ProgressReporter pr("PandaAnalyzer::Run",&iE,&nEvents,10);
	TimeReporter tr("PandaAnalyzer::Run",DEBUG);

	bool applyJSON = flags["applyJSON"];

	// EVENTLOOP --------------------------------------------------------------------------
	for (iE=nZero; iE!=nEvents; ++iE) {
		tr.Start();
		pr.Report();
		ResetBranches();
		tIn->GetEntry(iE);
		tr.TriggerEvent("GetEntry");

		// event info
		gt->mcWeight = (event->mcWeight>0) ? 1 : -1;
		gt->runNumber = event->runNumber;
		gt->lumiNumber = event->lumiNumber;
		gt->eventNumber = event->eventNumber;
		gt->npv = event->npv;
		gt->metFilter = (event->metfilters->at(0)) ? 1 : 0;
		if (!isData) 
			gt->sf_pu = hPUWeight->Eval(gt->npv);
		if (isData) {
			// check the json
			if (applyJSON && !PassGoodLumis(gt->runNumber,gt->lumiNumber)) 
				continue;

			// save triggers
			for (auto iT : metTriggers) {	 
				if (event->tiggers->at(iT)) {
					gt->trigger |= kMETTrig;
					break;
				}
			}
			for (auto iT : eleTriggers) {	 
				if (event->tiggers->at(iT)) {
					gt->trigger |= kSingleEleTrig;
					break;
				}
			}
			for (auto iT : phoTriggers) {	 
				if (event->tiggers->at(iT)) {
					gt->trigger |= kSinglePhoTrig;
					break;
				}
			}
		}

		tr.TriggerEvent("initialize");

		// default values for arrays
		for(unsigned int i=0;i<20;i++){
			gt->jetPt[i]=-99;
			gt->jetEta[i]=-99;
			gt->jetPhi[i]=-99;
			gt->jetE[i]=-99;
			gt->jetCSV[i]=-99;
			gt->jetIso[i]=-99;
			gt->jetQGL[i]=-99;
			if(i<2){
				gt->fj1sjPt[i]=-99;
				gt->fj1sjEta[i]=-99;
				gt->fj1sjPhi[i]=-99;
				gt->fj1sjM[i]=-99;
				gt->fj1sjCSV[i]=-99;
				gt->fj1sjQGL[i]=-99;
				gt->hbbjtidx[i]=-99;
			}
		}

		// met
		gt->pfmet = pfmet->pt;
		gt->pfmetphi = pfmet->phi;
		gt->calomet = pfmet->calo_pt;
		gt->puppimet = puppimet->pt;
		gt->puppimetphi = puppimet->phi;
		TLorentzVector vPFMET, vPuppiMET;
		vPFMET.SetPtEtaPhiM(gt->pfmet,0,gt->pfmetphi,0); 
		vPuppiMET.SetPtEtaPhiM(gt->puppimet,0,gt->puppimetphi,0); 
		TVector2 vMETNoMu; vMETNoMu.SetMagPhi(gt->pfmet,gt->pfmetphi); //	for trigger eff

		tr.TriggerEvent("met");

		//electrons
		std::vector<PObject*> looseLeps, tightLeps;
		for (PElectron *ele : *electrons) {
			float pt = ele->pt; float eta = ele->eta; float aeta = fabs(eta);
			if (pt<10 || aeta>2.5)
			//if (pt<10 || aeta>2.5 || (aeta>1.4442 && aeta<1.566))
				continue;
			if ((ele->id&PElectron::kVeto)==0)
				continue;
			// if (!ElectronIsolation(pt,eta,ele->iso,PElectron::kVeto))
			// 	continue;
			looseLeps.push_back(ele);
			gt->nLooseElectron++;
		} 

		// muons
		for (PMuon *mu : *muons) {
			float pt = mu->pt; float eta = mu->eta; float aeta = fabs(eta);
			if (pt<10 || aeta>2.4)
				continue;
			if ((mu->id&PMuon::kLoose)==0)
				continue;
			if (!MuonIsolation(pt,eta,mu->iso,PMuon::kLoose))
				continue;
			looseLeps.push_back(mu);
			gt->nLooseMuon++;
			TVector2 vMu; vMu.SetMagPhi(pt,mu->phi);
			vMETNoMu += vMu;
		} 
		gt->pfmetnomu = vMETNoMu.Mod();

		// now consider all leptons
		gt->nLooseLep = looseLeps.size();
		if (gt->nLooseLep>0) {
			int nToSort = TMath::Min(3,gt->nLooseLep);
			std::partial_sort(looseLeps.begin(),looseLeps.begin()+nToSort,looseLeps.end(),SortPObjects);
		}
		int lep_counter=1;
		for (PObject *lep : looseLeps) {
			if (lep_counter==1) {
				gt->looseLep1Pt = lep->pt;
				gt->looseLep1Eta = lep->eta;
				gt->looseLep1Phi = lep->phi;
			} else if (lep_counter==2) {
				gt->looseLep2Pt = lep->pt;
				gt->looseLep2Eta = lep->eta;
				gt->looseLep2Phi = lep->phi;
			} else {
				break;
			}
			// now specialize lepton types
			PMuon *mu = dynamic_cast<PMuon*>(lep);
			if (mu!=NULL) {
				bool isTight = ( (mu->id&PMuon::kTight)!=0 &&
													MuonIsolation(mu->pt,mu->eta,mu->iso,PMuon::kTight) &&
													mu->pt>20 && fabs(mu->eta)<2.4 );
				if (lep_counter==1) {
					gt->looseLep1PdgId = mu->q*-13;
					gt->looseLep1IsHLTSafe = 1;
					if (isTight) {
						gt->nTightMuon++;
						gt->looseLep1IsTight = 1;
						matchLeps.push_back(lep);
					}
				} else if (lep_counter==2) {
					gt->looseLep2PdgId = mu->q*-13;
					gt->looseLep2IsHLTSafe = 1;
					if (isTight) {
						gt->nTightMuon++;
						gt->looseLep2IsTight = 1;
					}
					if (isTight || gt->looseLep1IsTight)
						matchLeps.push_back(lep);
				}
			} else {
				PElectron *ele = dynamic_cast<PElectron*>(lep);
				bool isTight = ( (ele->id&PElectron::kTight)!=0 &&
													/*ElectronIsolation(ele->pt,ele->eta,ele->iso,PElectron::kTight) &&*/
													ele->pt>40 && fabs(ele->eta)<2.5 );
				if (lep_counter==1) {
					gt->looseLep1PdgId = ele->q*-11;
					gt->looseLep1IsHLTSafe = (ele->id&PElectron::kHLTPresel)!=0 ? 1 : 0;
					if (isTight) {
						gt->nTightElectron++;
						gt->looseLep1IsTight = 1;
						matchLeps.push_back(lep);
						matchEles.push_back(lep);
					}
				} else if (lep_counter==2) {
					gt->looseLep2PdgId = ele->q*-11;
					gt->looseLep2IsHLTSafe = (ele->id&PElectron::kHLTPresel)!=0 ? 1 : 0;
					if (isTight) {
						gt->nTightElectron++;
						gt->looseLep2IsTight = 1;
					}
					if (isTight || gt->looseLep1IsTight) {
						matchLeps.push_back(lep);
						matchEles.push_back(lep);
					}
				}
			}
			++lep_counter;
		}
		gt->nTightLep = gt->nTightElectron + gt->nTightMuon;
		if (gt->nLooseLep>0) {
			PObject *lep1 = looseLeps[0];
			gt->mT = MT(lep1->pt,lep1->phi,gt->pfmet,gt->pfmetphi);
		}
		if (gt->nLooseLep>1 && gt->looseLep1PdgId+gt->looseLep2PdgId==0) {
			TLorentzVector v1,v2;
			PObject *lep1=looseLeps[0], *lep2=looseLeps[1];
			v1.SetPtEtaPhiM(lep1->pt,lep1->eta,lep1->phi,lep1->m);
			v2.SetPtEtaPhiM(lep2->pt,lep2->eta,lep2->phi,lep2->m);
			gt->diLepMass = (v1+v2).M();
		} else {
			gt->diLepMass = -1;
		}

		tr.TriggerEvent("leptons"); 

		// photons
		std::vector<panda::PPhoton*> loosePhos;
		for (PPhoton *pho : *photons) {
			if ((pho->id&PPhoton::kLoose)==0 || (pho->id&PPhoton::kEleVeto)==0)
				continue;
			float pt = pho->pt;
			if (pt<1) continue;
			float eta = pho->eta, phi = pho->phi;
			if (pt<15 || fabs(eta)>2.5)
				continue;
			
			/*
			if (IsMatched(&matchEles,0.16,eta,phi))
				continue;
			*/
			
			loosePhos.push_back(pho);
			gt->nLoosePhoton++;
			if (gt->nLoosePhoton==1) {
				gt->loosePho1Pt = pt;
				gt->loosePho1Eta = eta;
				gt->loosePho1Phi = phi;
			}
			if ( (pho->id&PPhoton::kMedium)!=0 &&
						pt>175 /*&& fabs(eta)<1.4442*/ ) { // apply eta cut offline
				if (gt->nLoosePhoton==1)
					gt->loosePho1IsTight=1;
				gt->nTightPhoton++;
				matchPhos.push_back(pho);
			}
		}

		if (isData && gt->nLoosePhoton>0) {
			if (gt->loosePho1Pt>=175 && gt->loosePho1Pt<200)
				gt->sf_phoPurity = 0.04802;
			else if (gt->loosePho1Pt>=200 && gt->loosePho1Pt<250)
				gt->sf_phoPurity = 0.04241;
			else if (gt->loosePho1Pt>=250 && gt->loosePho1Pt<300)
				gt->sf_phoPurity = 0.03641;
			else if (gt->loosePho1Pt>=300 && gt->loosePho1Pt<350)
				gt->sf_phoPurity = 0.0333;
			else if (gt->loosePho1Pt>=350)
				gt->sf_phoPurity = 0.02544;
		}

		tr.TriggerEvent("photons");

		// trigger efficiencies
		gt->sf_eleTrig=1; gt->sf_metTrig=1; gt->sf_phoTrig=1;
		if (!isData) {
			gt->sf_metTrig = hMetTrig->Eval(gt->pfmetnomu);

			if (gt->nLooseElectron>0 && abs(gt->looseLep1PdgId)==1
					&& gt->looseLep1IsTight==1 && gt->nLooseMuon==0) {
				float eff1=0, eff2=0;
				if (gt->looseLep1Pt<100) {
					eff1 = hEleTrigLow->Eval(gt->looseLep1Eta,gt->looseLep1Pt);
				} else {
					if (fabs(gt->looseLep1Eta)<1.4442) {
						eff1 = hEleTrigB->Eval(gt->looseLep1Pt);
					}
					if (1.5660<fabs(gt->looseLep1Eta) && fabs(gt->looseLep1Eta)<2.5) {
						eff1 = hEleTrigE->Eval(gt->looseLep1Pt);
					}
				}
				if (gt->nLooseElectron>1 && fabs(gt->looseLep2PdgId)==11) {
					if (gt->looseLep2Pt<100) {
						eff2 = hEleTrigLow->Eval(gt->looseLep2Eta,gt->looseLep2Pt);
					} else {
						if (fabs(gt->looseLep2Eta)<1.4442) {
							eff2 = hEleTrigB->Eval(gt->looseLep2Pt);
						}
						if (1.5660<fabs(gt->looseLep2Eta) && fabs(gt->looseLep2Eta)<2.5) {
							eff2 = hEleTrigE->Eval(gt->looseLep2Pt);
						}
					}
				}
				gt->sf_eleTrig = 1 - (1-eff1)*(1-eff2);
			} // done with ele trig SF

			if (gt->nLoosePhoton>0 && gt->loosePho1IsTight)
				gt->sf_phoTrig = hPhoTrig->Eval(gt->loosePho1Pt);
		} 

		tr.TriggerEvent("triggers");

		// recoil!
		TLorentzVector vObj1, vObj2;
		TLorentzVector vUW, vUZ, vUA;
		TLorentzVector vpfUW, vpfUZ, vpfUA;
		if (gt->nLooseLep>0) {
			PObject *lep1 = looseLeps.at(0);
			vObj1.SetPtEtaPhiM(lep1->pt,lep1->eta,lep1->phi,lep1->m);

			// one lep => W
			vUW = vPuppiMET+vObj1; gt->UWmag=vUW.Pt(); gt->UWphi=vUW.Phi();
			vpfUW = vPFMET+vObj1; gt->pfUWmag=vpfUW.Pt(); gt->pfUWphi=vpfUW.Phi();

			if (gt->nLooseLep>1 && gt->looseLep1PdgId+gt->looseLep2PdgId==0) {
				// two OS lep => Z
				PObject *lep2 = looseLeps.at(1);
				vObj2.SetPtEtaPhiM(lep2->pt,lep2->eta,lep2->phi,lep2->m);

				vUZ=vUW+vObj2; gt->UZmag=vUZ.Pt(); gt->UZphi=vUZ.Phi();
				vpfUZ=vpfUW+vObj2; gt->pfUZmag=vpfUZ.Pt(); gt->pfUZphi=vpfUZ.Phi();
			}
		}
		if (gt->nLoosePhoton>0) {
			PPhoton *pho = loosePhos.at(0);
			vObj1.SetPtEtaPhiM(pho->pt,pho->eta,pho->phi,pho->m);

			vUA=vPuppiMET+vObj1; gt->UAmag=vUA.Pt(); gt->UAphi=vUA.Phi();
			vpfUA=vPFMET+vObj1; gt->pfUAmag=vpfUA.Pt(); gt->pfUAphi=vpfUA.Phi();
		}

		tr.TriggerEvent("recoils");

		PFatJet *fj1=0;
		gt->nFatjet=0;
		if (flags["fatjet"]) {
			int fatjet_counter=-1;
			for (PFatJet *fj : *fatjets) {
				++fatjet_counter;
				float pt = fj->pt;
				float rawpt = fj->rawPt;
				float eta = fj->eta;
				float mass = fj->m;
				float ptcut = 200;
				if (flags["monohiggs"]) 
					ptcut = 200;
				if (pt<ptcut || fabs(eta)>2.4 || (fj->id&PFatJet::kMonojet)==0)
					continue;
	
				float phi = fj->phi;
				if (IsMatched(&matchLeps,2.25,eta,phi) || IsMatched(&matchPhos,2.25,eta,phi)) {
					continue;
				}
				gt->nFatjet++;
				if (gt->nFatjet==1) {
					fj1 = fj;
					if (fatjet_counter==0)
						gt->fj1IsClean = 1;
					else
						gt->fj1IsClean = 0;
					gt->fj1Pt = pt;
					gt->fj1Eta = eta;
					gt->fj1Phi = phi;
					gt->fj1M = mass;
					gt->fj1MSD = fj->mSD;
					gt->fj1RawPt = rawpt;

					// do a bit of jet energy scaling 
					JetCorrectionUncertainty *uncReader=0;
					if (isData) {
						TString thisEra = eras.getEra(gt->runNumber);
						for (auto &iter : ak8UncReader) {
							if (! iter.first.Contains("data"))
								continue;
							if (iter.first.Contains(thisEra)) {
								uncReader = iter.second;
								break;
							}
						}
					} else {
						uncReader = ak8UncReader["MC"];
					}
					uncReader->setJetEta(eta); uncReader->setJetPt(pt);
					double scaleUnc = uncReader->getUncertainty(true);
					gt->fj1PtScaleUp    = gt->fj1Pt  * (1 + 2*scaleUnc);
					gt->fj1PtScaleDown  = gt->fj1Pt  * (1 - 2*scaleUnc);
					gt->fj1MSDScaleUp   = gt->fj1MSD * (1 + 2*scaleUnc);
					gt->fj1MSDScaleDown = gt->fj1MSD * (1 - 2*scaleUnc);

					// do some jet energy smearing
					if (isData) {
						gt->fj1PtSmeared = gt->fj1Pt;
						gt->fj1PtSmearedUp = gt->fj1Pt;
						gt->fj1PtSmearedDown = gt->fj1Pt;
						gt->fj1MSDSmeared = gt->fj1MSD;
						gt->fj1MSDSmearedUp = gt->fj1MSD;
						gt->fj1MSDSmearedDown = gt->fj1MSD;
					} else {
						double smear=1, smearUp=1, smearDown=1;
					  ak8JERReader->getStochasticSmear(pt,eta,15,smear,smearUp,smearDown);

						gt->fj1PtSmeared = smear*gt->fj1Pt;
						gt->fj1PtSmearedUp = smearUp*gt->fj1Pt;
						gt->fj1PtSmearedDown = smearDown*gt->fj1Pt;

						gt->fj1MSDSmeared = smear*gt->fj1MSD;
						gt->fj1MSDSmearedUp = smearUp*gt->fj1MSD;
						gt->fj1MSDSmearedDown = smearDown*gt->fj1MSD;
					}

					// mSD correction
					if (flags["monohiggs"]) {
						float corrweight=1.;
						corrweight = getMSDcorr(pt,eta);
						gt->fj1MSD_corr = corrweight*gt->fj1MSD;
					}

					// now we do substructure
					gt->fj1Tau32 = clean(fj->tau3/fj->tau2);
					gt->fj1Tau32SD = clean(fj->tau3SD/fj->tau2SD);
					gt->fj1Tau21 = clean(fj->tau2/fj->tau1);
					gt->fj1Tau21SD = clean(fj->tau2SD/fj->tau1SD);

					for (unsigned int iB=0; iB!=betas.size(); ++iB) {
						float beta = betas.at(iB);
						for (auto N : Ns) {
							for (auto order : orders) {
								if (gt->fj1IsClean || true)
									gt->fj1ECFNs[makeECFString(order,N,beta)] = fj->get_ecf(order,N,iB); 
								else
									gt->fj1ECFNs[makeECFString(order,N,beta)] = -1; 
							}
						}
					} //loop over betas
					gt->fj1HTTMass = fj->htt_mass;
					gt->fj1HTTFRec = fj->htt_frec;

					VJet *subjets = fj->subjets;

					std::sort(subjets->begin(),subjets->end(),SortPJetByCSV);
					gt->fj1MaxCSV = subjets->at(0)->csv; 
					gt->fj1MinCSV = subjets->back()->csv; 

					if (flags["monohiggs"]) {
						for (unsigned int iSJ=0; iSJ!=fj->subjets->size(); ++iSJ) {
							PJet *subjet = fj1->subjets->at(iSJ);
							gt->fj1sjPt[iSJ]=subjet->pt;
							gt->fj1sjEta[iSJ]=subjet->eta;
							gt->fj1sjPhi[iSJ]=subjet->phi;
							gt->fj1sjM[iSJ]=subjet->m;
							gt->fj1sjCSV[iSJ]=subjet->csv;
							gt->fj1sjQGL[iSJ]=subjet->qgl;
						}
					}
				}
			}
		}

		tr.TriggerEvent("fatjet");

		// first identify interesting jets
		vector<PJet*> cleanedJets, isoJets, btaggedJets, centralJets;
		vector<int> btagindices;
		TLorentzVector vJet;
		PJet *jet1=0, *jet2=0;
		gt->dphipuppimet=999; gt->dphipfmet=999;
		gt->dphiUW=999; gt->dphipfUW=999;
		gt->dphiUZ=999; gt->dphipfUZ=999;
		gt->dphiUA=999; gt->dphipfUA=999;
		for (PJet *jet : *jets) {
			if (jet->pt<30 || abs(jet->eta)>4.5) 
				continue;
			if ((jet->id&PJet::kLoose)==0)
				continue; // apply loose selection to all jets
			if (IsMatched(&matchLeps,0.16,jet->eta,jet->phi) ||
					IsMatched(&matchPhos,0.16,jet->eta,jet->phi))
				continue;

			cleanedJets.push_back(jet);
			float csv = (fabs(jet->eta)<2.5) ? jet->csv : -1;
			if (fabs(jet->eta)<2.4) {
				centralJets.push_back(jet);
				if (centralJets.size()==1) {
					jet1 = jet;
					gt->jet1Pt = jet->pt;
					gt->jet1Eta = jet->eta;
					gt->jet1Phi = jet->phi;
					gt->jet1CSV = csv; 
					gt->jet1IsTight = ((jet->id&PJet::kMonojet)==0) ? 0 : 1;
				} else if (centralJets.size()==2) {
					jet2 = jet;
					gt->jet2Pt = jet->pt;
					gt->jet2Eta = jet->eta;
					gt->jet2Phi = jet->phi;
					gt->jet2CSV = csv; 
				} 
			}

			if (flags["monohiggs"]) {
				gt->jetPt[cleanedJets.size()-1]=jet->pt;
				gt->jetEta[cleanedJets.size()-1]=jet->eta;
				gt->jetPhi[cleanedJets.size()-1]=jet->phi;
				gt->jetE[cleanedJets.size()-1]=jet->m;
				gt->jetCSV[cleanedJets.size()-1]=csv;
				gt->jetQGL[cleanedJets.size()-1]=jet->qgl;
			}

			// compute dphi wrt mets
			if (cleanedJets.size()<5) {
				vJet.SetPtEtaPhiM(jet->pt,jet->eta,jet->phi,jet->m);
				gt->dphipuppimet = std::min(fabs(vJet.DeltaPhi(vPuppiMET)),(double)gt->dphipuppimet);
				gt->dphipfmet = std::min(fabs(vJet.DeltaPhi(vPFMET)),(double)gt->dphipfmet);
				gt->dphiUA = std::min(fabs(vJet.DeltaPhi(vUA)),(double)gt->dphiUA);
				gt->dphiUW = std::min(fabs(vJet.DeltaPhi(vUW)),(double)gt->dphiUW);
				gt->dphiUZ = std::min(fabs(vJet.DeltaPhi(vUZ)),(double)gt->dphiUZ);
				gt->dphipfUA = std::min(fabs(vJet.DeltaPhi(vpfUA)),(double)gt->dphipfUA);
				gt->dphipfUW = std::min(fabs(vJet.DeltaPhi(vpfUW)),(double)gt->dphipfUW);
				gt->dphipfUZ = std::min(fabs(vJet.DeltaPhi(vpfUZ)),(double)gt->dphipfUZ);
			}
			// btags
			if (csv>0.8484) { 
				++gt->jetNBtags;
				if (flags["monohiggs"]) {
					btaggedJets.push_back(jet);
					btagindices.push_back(cleanedJets.size()-1);
				}
			}

			bool isIsoJet = false;
			if (gt->nFatjet==0) {
				isIsoJet = true;
			} else if (fabs(jet->eta)<2.5 
					       && DeltaR2(gt->fj1Eta,gt->fj1Phi,jet->eta,jet->phi)>2.25) {
				isIsoJet = true;
			
			}

			if (isIsoJet) {
				isoJets.push_back(jet);
				if (csv>0.5426) 
					++gt->isojetNBtags;
				if (isoJets.size()==1) {
					gt->isojet1Pt = jet->pt;
					gt->isojet1CSV = jet->csv;
				} else if (isoJets.size()==2) {
					gt->isojet2Pt = jet->pt;
					gt->isojet2CSV = jet->csv;
				}
				if (flags["monohiggs"]) 
					gt->jetIso[cleanedJets.size()-1]=1;
			} else {
				if (flags["monohiggs"]) 
					gt->jetIso[cleanedJets.size()-1]=0;
			}

		} // VJet loop

		gt->nIsoJet = isoJets.size();
		gt->nJet = cleanedJets.size();
		if (gt->nJet>1 && flags["monojet"]) {
			gt->jet12DEta = fabs(jet1->eta-jet2->eta);
			TLorentzVector vj1, vj2; 
			vj1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->m);
			vj2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->m);
			gt->jet12Mass = (vj1+vj2).M();
		}

		tr.TriggerEvent("jets");


		if (flags["monohiggs"]){
			// Higgs reconstrcution for resolved analysis - highest pt pair of b jets
			float tmp_hbbpt=-99;
			float tmp_hbbeta=-99;
			float tmp_hbbphi=-99;
			float tmp_hbbm=-99;
			int tmp_hbbjtidx1=-1;
			int tmp_hbbjtidx2=-1;
			for (unsigned int i = 0;i<btaggedJets.size();i++){
				PJet *jet_1 = btaggedJets.at(i);
				TLorentzVector hbbdaughter1;
				hbbdaughter1.SetPtEtaPhiM(jet_1->pt,jet_1->eta,jet_1->phi,jet_1->m);
				for (unsigned int j = i+1;j<btaggedJets.size();j++){
					PJet *jet_2 = btaggedJets.at(j);
					TLorentzVector hbbdaughter2;
					hbbdaughter2.SetPtEtaPhiM(jet_2->pt,jet_2->eta,jet_2->phi,jet_2->m);
					TLorentzVector hbbsystem = hbbdaughter1 + hbbdaughter2;
					if (hbbsystem.Pt()>tmp_hbbpt){
						tmp_hbbpt = hbbsystem.Pt();
						tmp_hbbeta = hbbsystem.Eta();
						tmp_hbbphi = hbbsystem.Phi();
						tmp_hbbm = hbbsystem.M();
						tmp_hbbjtidx1 = btagindices.at(i);
						tmp_hbbjtidx2 = btagindices.at(j);
					}
				}
			}
			gt->hbbpt = tmp_hbbpt;
			gt->hbbeta = tmp_hbbeta;
			gt->hbbphi = tmp_hbbphi;
			gt->hbbm = tmp_hbbm;
			gt->hbbjtidx[0] = tmp_hbbjtidx1;
			gt->hbbjtidx[1] = tmp_hbbjtidx2;

			tr.TriggerEvent("monohiggs");
		}

		for (PTau *tau : *taus) {
			if ((tau->id&PTau::kDecayModeFinding)==0 ||
					(tau->id&PTau::kDecayModeFindingNewDMs)==0) 
				continue;
			if ((tau->id&PTau::kVLooseIsolationMVArun2v1DBnewDMwLT)==0)
				continue;
//			if (tau->isoDeltaBetaCorr>5)
//				continue; 
			if (tau->pt<18 || fabs(tau->eta)>2.3)
				continue;
			if (IsMatched(&matchLeps,0.16,tau->eta,tau->phi))
				continue;
			gt->nTau++;
		} 

		tr.TriggerEvent("taus");

		if (!PassPreselection())
			continue;

		tr.TriggerEvent("presel");

		// identify interesting gen particles for fatjet matching
		unsigned int pdgidTarget=0;
		if (!isData && processType>=kTT) {
			switch(processType) {
				case kTop:
				case kTT:
					pdgidTarget=6;
					break;
				case kV:
					pdgidTarget=24;
					break;
				case kH:
					pdgidTarget=25;
					break;
				default:
					// processType>=kTT means we should never get here
					PError("PandaAnalyzer::Run","Reached an unknown process type");
			}

			std::vector<int> targets;

			int nGen = genparts->size();
			for (int iG=0; iG!=nGen; ++iG) {
				PGenParticle *part = genparts->at(iG);
				int pdgid = part->pdgid;
				unsigned int abspdgid = abs(pdgid);
				if (abspdgid == pdgidTarget)
					targets.push_back(iG);
			} //looking for targets

			for (int iG : targets) {
				PGenParticle *part = genparts->at(iG);

				// check there is no further copy:
				bool isLastCopy=true;
				for (int jG : targets) {
					if (genparts->at(jG)->parent==iG) {
						isLastCopy=false;
						break;
					}
				}
				if (!isLastCopy)
					continue;
				
				// (a) check it is a hadronic decay and if so, (b) calculate the size
				if (processType==kTop||processType==kTT) {
					
					// first look for a W whose parent is the top at iG, or a W further down the chain
					int iW=-1;
					for (int jG=0; jG!=nGen; ++jG) {
						PGenParticle *partW = genparts->at(jG);
						if (TMath::Abs(partW->pdgid)==24 && partW->pdgid*part->pdgid>0) {
							// it's a W and has the same sign as the top
							if (iW<0 && partW->parent==iG) {
								iW=jG;
							} else if (iW>=0 && partW->parent==iW) {
								iW=jG;
							}
						}
					} // looking for W
					if (iW<0) {// ???
						continue;
					}
					PGenParticle *partW = genparts->at(iW);

					// now look for b or W->qq
					int iB=-1, iQ1=-1, iQ2=-1;
					double size=0, sizeW=0;
					for (int jG=0; jG!=nGen; ++jG) {
						PGenParticle *partQ = genparts->at(jG);
						int pdgidQ = partQ->pdgid;
						unsigned int abspdgidQ = TMath::Abs(pdgidQ);
						if (abspdgidQ>5)
							continue;
						if (abspdgidQ==5 && iB<0 && partQ->parent==iG) {
							// only keep first copy
							iB = jG;
							size = TMath::Max(DeltaR2(part->eta,part->phi,partQ->eta,partQ->phi),size);
						} else if (abspdgidQ<5 && partQ->parent==iW) {
							if (iQ1<0) {
								iQ1 = jG;
								size = TMath::Max(DeltaR2(part->eta,part->phi,partQ->eta,partQ->phi),size);
								sizeW = TMath::Max(DeltaR2(partW->eta,partW->phi,partQ->eta,partQ->phi),sizeW);
							} else if (iQ2<0) {
								iQ2 = jG;
								size = TMath::Max(DeltaR2(part->eta,part->phi,partQ->eta,partQ->phi),size);
								sizeW = TMath::Max(DeltaR2(partW->eta,partW->phi,partQ->eta,partQ->phi),sizeW);
							}
						}
						if (iB>=0 && iQ1>=0 && iQ2>=0)
							break;
					} // looking for quarks


					bool isHadronic = (iB>=0 && iQ1>=0 && iQ2>=0); // all 3 quarks were found
					if (isHadronic)
						genObjects[part] = size;

					bool isHadronicW = (iQ1>=0 && iQ2>=0);
					if (isHadronicW)
						genObjects[partW] = sizeW;

				} else { // these are W,Z,H - 2 prong decays

					int iQ1=-1, iQ2=-1;
					double size=0;
					for (int jG=0; jG!=nGen; ++jG) {
						PGenParticle *partQ = genparts->at(jG);
						int pdgidQ = partQ->pdgid;
						unsigned int abspdgidQ = TMath::Abs(pdgidQ);
						if (abspdgidQ>5)
							continue; 
						if (partQ->parent==iG) {
							if (iQ1<0) {
								iQ1=jG;
								size = TMath::Max(DeltaR2(part->eta,part->phi,partQ->eta,partQ->phi),size);
							} else if (iQ2<0) {
								iQ2=jG;
								size = TMath::Max(DeltaR2(part->eta,part->phi,partQ->eta,partQ->phi),size);
							}
						}
						if (iQ1>=0 && iQ2>=0)
							break;
					} // looking for quarks

					bool isHadronic = (iQ1>=0 && iQ2>=0); // both quarks were found

					// add to coll0ection
					if (isHadronic)
						genObjects[part] = size;
				}

			} // loop over targets
		} // process is interesting

		tr.TriggerEvent("gen matching");

		// do gen matching now that presel is passed
		if (!isData && gt->nFatjet>0) {
			// first see if jet is matched
			PGenParticle *matched = MatchToGen(fj1->eta,fj1->phi,1.5,pdgidTarget);
			if (matched!=NULL) {
				gt->fj1IsMatched = 1;
				gt->fj1GenPt = matched->pt;
				gt->fj1GenSize = genObjects[matched];
			} else {
				gt->fj1IsMatched = 0;
			}
			if (pdgidTarget==6) { // matched to top; try for W
				PGenParticle *matchedW = MatchToGen(fj1->eta,fj1->phi,1.5,24);
				if (matchedW!=NULL) {
					gt->fj1IsWMatched = 1;
					gt->fj1GenWPt = matchedW->pt;
					gt->fj1GenWSize = genObjects[matchedW];
				} else {
					gt->fj1IsWMatched = 0;
				}
			}
			
			// now get the highest pT gen particle inside the jet cone
			for (PGenParticle *gen : *genparts) {
				float pt = gen->pt;
				int pdgid = gen->pdgid;
				if (pt>(gt->fj1HighestPtGenPt) 
						&& DeltaR2(gen->eta,gen->phi,fj1->eta,fj1->phi)<2.25) {
					gt->fj1HighestPtGenPt = pt;
					gt->fj1HighestPtGen = pdgid;
				}
			} 

			// now get the subjet btag SFs
			vector<btagcand> sj_btagcands;
			vector<double> sj_sf_cent, sj_sf_bUp, sj_sf_bDown, sj_sf_mUp, sj_sf_mDown;
			unsigned int nSJ = fj1->subjets->size();
			for (unsigned int iSJ=0; iSJ!=nSJ; ++iSJ) {
				PJet *subjet = fj1->subjets->at(iSJ);
				int flavor=0;
				for (PGenParticle *gen : *genparts) {
					int apdgid = abs(gen->pdgid);
					if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
						continue;
					double dr2 = DeltaR2(subjet->eta,subjet->phi,gen->eta,gen->phi);
					if (dr2<0.09) {
						if (apdgid==4 || apdgid==5) {
							flavor=apdgid;
							break;
						} else {
							flavor=0;
						}
					}
				} // finding the subjet flavor
	
				float pt = subjet->pt;
				float btagUncFactor = 1;
				float eta = subjet->eta;
				double eff(1),sf(1),sfUp(1),sfDown(1);
				unsigned int binpt = btagpt.bin(pt);
				unsigned int bineta = btageta.bin(fabs(eta));
				if (flavor==5) {
					eff = beff[bineta][binpt];
				} else if (flavor==4) {
					eff = ceff[bineta][binpt];
				} else {
					eff = lfeff[bineta][binpt];
				}
				calcBJetSFs("sj_L",flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
				sj_btagcands.push_back(btagcand(iSJ,flavor,eff,sf,sfUp,sfDown));
				sj_sf_cent.push_back(sf);
				if (flavor>0) {
					sj_sf_bUp.push_back(sfUp); sj_sf_bDown.push_back(sfDown);
					sj_sf_mUp.push_back(sf); sj_sf_mDown.push_back(sf);
				} else {
					sj_sf_bUp.push_back(sf); sj_sf_bDown.push_back(sf);
					sj_sf_mUp.push_back(sfUp); sj_sf_mDown.push_back(sfDown);
				}

				// evaluate the CSV weight
				if (subjet->csv>0) {
					if (flavor==0) {
						gt->sf_sjcsvWeightM *= hCSVLF->Eval(subjet->csv);
					} else {
						gt->sf_sjcsvWeightB *= hCSVHF->Eval(subjet->csv);
					}
				}
	
			} // loop over subjets
			gt->sf_sjcsvWeightMUp = 1.05 * gt->sf_sjcsvWeightM;
			gt->sf_sjcsvWeightMDown = 0.95 * gt->sf_sjcsvWeightM;
			gt->sf_sjcsvWeightBUp = 1.05 * gt->sf_sjcsvWeightB;
			gt->sf_sjcsvWeightBDown = 0.95 * gt->sf_sjcsvWeightB;

			EvalBtagSF(sj_btagcands,sj_sf_cent,
									gt->sf_sjbtag0,gt->sf_sjbtag1,gt->sf_sjbtag2);
			EvalBtagSF(sj_btagcands,sj_sf_bUp,
									gt->sf_sjbtag0BUp,gt->sf_sjbtag1BUp,gt->sf_sjbtag2BUp);
			EvalBtagSF(sj_btagcands,sj_sf_bDown,
									gt->sf_sjbtag0BDown,gt->sf_sjbtag1BDown,gt->sf_sjbtag2BDown);
			EvalBtagSF(sj_btagcands,sj_sf_mUp,
									gt->sf_sjbtag0MUp,gt->sf_sjbtag1MUp,gt->sf_sjbtag2MUp);
			EvalBtagSF(sj_btagcands,sj_sf_mDown,
									gt->sf_sjbtag0MDown,gt->sf_sjbtag1MDown,gt->sf_sjbtag2MDown);

		}

		tr.TriggerEvent("fatjet gen-matching");

		if (!isData) {
			// now get the jet btag SFs
			vector<btagcand> btagcands;
			vector<btagcand> btagcands_alt;
			vector<double> sf_cent, sf_bUp, sf_bDown, sf_mUp, sf_mDown;
			vector<double> sf_cent_alt, sf_bUp_alt, sf_bDown_alt, sf_mUp_alt, sf_mDown_alt;

			unsigned int nJ = centralJets.size();
			for (unsigned int iJ=0; iJ!=nJ; ++iJ) {
				PJet *jet = centralJets.at(iJ);
				bool isIsoJet=false;
				if (std::find(isoJets.begin(), isoJets.end(), jet) != isoJets.end())
					isIsoJet = true;
				int flavor=0;
				float genpt=0;
				for (PGenParticle *gen : *genparts) {
					int apdgid = abs(gen->pdgid);
					if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
						continue;
					double dr2 = DeltaR2(jet->eta,jet->phi,gen->eta,gen->phi);
					if (dr2<0.09) {
						genpt = gen->pt;
						if (apdgid==4 || apdgid==5) {
							flavor=apdgid;
							break;
						} else {
							flavor=0;
						}
					}
				} // finding the jet flavor
				float pt = jet->pt;
				float btagUncFactor = 1;
				float eta = jet->eta;
				double eff(1),sf(1),sfUp(1),sfDown(1);
				unsigned int binpt = btagpt.bin(pt);
				unsigned int bineta = btageta.bin(fabs(eta));
				if (flavor==5) 
					eff = beff[bineta][binpt];
				else if (flavor==4) 
					eff = ceff[bineta][binpt];
				else 
					eff = lfeff[bineta][binpt];
				if (jet==centralJets.at(0)) {
					gt->jet1Flav = flavor;
					gt->jet1GenPt = genpt;
				} else if (jet==centralJets.at(1)) {
					gt->jet2Flav = flavor;
					gt->jet2GenPt = genpt;
				}
				if (isIsoJet) {
					if (jet==isoJets.at(0)) 
						gt->isojet1Flav = flavor;
					else if (jet==isoJets.at(1))
						gt->isojet2Flav = flavor;

					calcBJetSFs("jet_L",flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
					btagcands.push_back(btagcand(iJ,flavor,eff,sf,sfUp,sfDown));
					sf_cent.push_back(sf);

					if (flavor>0) {
						sf_bUp.push_back(sfUp); sf_bDown.push_back(sfDown);
						sf_mUp.push_back(sf); sf_mDown.push_back(sf);
					} else {
						sf_bUp.push_back(sf); sf_bDown.push_back(sf);
						sf_mUp.push_back(sfUp); sf_mDown.push_back(sfDown);
					}

					// evaluate the CSV weight
					if (jet->csv>0) {
						if (flavor==0) {
							gt->sf_csvWeightM *= hCSVLF->Eval(jet->csv);
						} else {
							gt->sf_csvWeightB *= hCSVHF->Eval(jet->csv);
						}
					}

				}

				// alternate stuff for inclusive jet collection (also different b tagging WP)
				double sf_alt(1),sfUp_alt(1),sfDown_alt(1); 
				calcBJetSFs("jet_M",flavor,eta,pt,eff,btagUncFactor,sf_alt,sfUp_alt,sfDown_alt);
				btagcands_alt.push_back(btagcand(iJ,flavor,eff,sf_alt,sfUp_alt,sfDown_alt));
				sf_cent_alt.push_back(sf_alt);
				if (flavor>0) {
					sf_bUp_alt.push_back(sfUp_alt); sf_bDown_alt.push_back(sfDown_alt);
					sf_mUp_alt.push_back(sf_alt); sf_mDown_alt.push_back(sf_alt);
				} else {
					sf_bUp_alt.push_back(sf_alt); sf_bDown_alt.push_back(sf_alt);
					sf_mUp_alt.push_back(sfUp_alt); sf_mDown_alt.push_back(sfDown_alt);
				}
			} // loop over jets

			gt->sf_csvWeightMUp = 1.05 * gt->sf_csvWeightM;
			gt->sf_csvWeightMDown = 0.95 * gt->sf_csvWeightM;
			gt->sf_csvWeightBUp = 1.05 * gt->sf_csvWeightB;
			gt->sf_csvWeightBDown = 0.95 * gt->sf_csvWeightB;

			EvalBtagSF(btagcands,sf_cent,
									gt->sf_btag0,gt->sf_btag1,gt->sf_btag2);
			EvalBtagSF(btagcands,sf_bUp,
									gt->sf_btag0BUp,gt->sf_btag1BUp,gt->sf_btag2BUp);
			EvalBtagSF(btagcands,sf_bDown,
									gt->sf_btag0BDown,gt->sf_btag1BDown,gt->sf_btag2BDown);
			EvalBtagSF(btagcands,sf_mUp,
									gt->sf_btag0MUp,gt->sf_btag1MUp,gt->sf_btag2MUp);
			EvalBtagSF(btagcands,sf_mDown,
									gt->sf_btag0MDown,gt->sf_btag1MDown,gt->sf_btag2MDown);

			EvalBtagSF(btagcands_alt,sf_cent_alt,
								 gt->sf_btag0_alt,gt->sf_btag1_alt,gt->sf_btag2_alt);
			EvalBtagSF(btagcands_alt,sf_bUp_alt,
								 gt->sf_btag0BUp_alt,gt->sf_btag1BUp_alt,gt->sf_btag2BUp_alt);
			EvalBtagSF(btagcands_alt,sf_bDown_alt,
								 gt->sf_btag0BDown_alt,gt->sf_btag1BDown_alt,gt->sf_btag2BDown_alt);
			EvalBtagSF(btagcands_alt,sf_mUp_alt,
								 gt->sf_btag0MUp_alt,gt->sf_btag1MUp_alt,gt->sf_btag2MUp_alt);
			EvalBtagSF(btagcands_alt,sf_mDown_alt,
								 gt->sf_btag0MDown_alt,gt->sf_btag1MDown_alt,gt->sf_btag2MDown_alt);
		}

		tr.TriggerEvent("ak4 gen-matching");

		// ttbar pT weight and NLO-NNLO kfactor
		// both correct for the same effect, but one comes from data
		// and the other from calculation
		gt->sf_tt = 1; gt->sf_tt_ext = 1; gt->sf_tt_bound = 1;
		gt->sf_tt8TeV = 1; gt->sf_tt8TeV_ext = 1; gt->sf_tt8TeV_bound = 1;
		gt->sf_qcdTT = 0;
		if (!isData && processType==kTT) {
			gt->genWPlusPt = -1; gt->genWMinusPt = -1;
			for (auto *gen : *genparts) {
				if (abs(gen->pdgid)!=24)
					continue;
				if (flags["firstGen"]) {
					if (gen->parent>=0 && genparts->at(gen->parent)->pdgid==gen->pdgid)
						continue; // must be first copy
				}
				if (gen->pdgid>0) {
					gt->genWPlusPt = gen->pt;
					gt->genWPlusEta = gen->eta;
				} else {
					gt->genWMinusPt = gen->pt;
					gt->genWMinusEta = gen->eta;
				}
				if (flags["firstGen"]) {
					if (gt->genWPlusPt>0 && gt->genWMinusPt>0)
						break;
				}
			}
			TLorentzVector vT,vTbar;
			float pt_t=0, pt_tbar=0;
			for (auto *gen : *genparts) {
				if (abs(gen->pdgid)!=6)
					continue;
				if (flags["firstGen"]) {
					if (gen->parent>=0 && genparts->at(gen->parent)->pdgid==gen->pdgid)
						continue; // must be first copy
				}
				if (gen->pdgid>0) {
					pt_t = gen->pt;
					gt->genTopPt = gen->pt;
					gt->genTopEta = gen->eta;
					vT.SetPtEtaPhiM(gen->pt,gen->eta,gen->phi,gen->m);
				} else {
					pt_tbar = gen->pt;
					gt->genAntiTopPt = gen->pt;
					gt->genAntiTopEta = gen->eta;
					vTbar.SetPtEtaPhiM(gen->pt,gen->eta,gen->phi,gen->m);
				}
				if (flags["firstGen"]) {
					if (pt_t>0 && pt_tbar>0)
						break;
				}
			}
			if (pt_t>0 && pt_tbar>0) {
				TLorentzVector vTT = vT+vTbar;
				gt->genTTPt = vTT.Pt(); gt->genTTEta = vTT.Eta();
				gt->sf_tt8TeV = TMath::Sqrt( 
                          TMath::Exp(0.156-0.00137*TMath::Min((float)400.,pt_t)) * 
                          TMath::Exp(0.156-0.00137*TMath::Min((float)400.,pt_tbar)) 
                          );
				gt->sf_tt = TMath::Sqrt( 
                          TMath::Exp(0.0615-0.0005*TMath::Min((float)400.,pt_t)) * 
                          TMath::Exp(0.0615-0.0005*TMath::Min((float)400.,pt_tbar)) 
                          );
				gt->sf_tt8TeV_ext = TMath::Sqrt( 
                          TMath::Exp(0.156-0.00137*pt_t) * 
                          TMath::Exp(0.156-0.00137*pt_tbar) 
                          );
				gt->sf_tt_ext = TMath::Sqrt( 
                          TMath::Exp(0.0615-0.0005*pt_t) * 
                          TMath::Exp(0.0615-0.0005*pt_tbar) 
                          );
				gt->sf_tt8TeV_bound = TMath::Sqrt( 
                          ((pt_t>400) ? 1 : TMath::Exp(0.156-0.00137*pt_t)) * 
                          ((pt_tbar>400) ? 1 : TMath::Exp(0.156-0.00137*pt_tbar))
                          );
				gt->sf_tt_bound = TMath::Sqrt( 
                          ((pt_t>400) ? 1 : TMath::Exp(0.0615-0.0005*pt_t)) * 
                          ((pt_tbar>400) ? 1 : TMath::Exp(0.0615-0.0005*pt_tbar))
                          );
			}

			if (pt_t>0)
				gt->sf_qcdTT *= TTNLOToNNLO(pt_t);
			if (pt_tbar>0) 
				gt->sf_qcdTT *= TTNLOToNNLO(pt_tbar);

			gt->sf_qcdTT = TMath::Sqrt(gt->sf_qcdTT);

			tr.TriggerEvent("tt SFs");
		} 

		// derive ewk/qcd weights
		gt->sf_qcdV=1; gt->sf_ewkV=1;
		if (!isData) {
			bool found = processType!=kA && processType!=kZ && processType!=kW;
			int target=24;
			if (processType==kZ) target=23;
			if (processType==kA) target=22;
			for (PGenParticle *gen : *genparts) {
				if (found) break;
				int apdgid = abs(gen->pdgid);
				if (apdgid==target)	{
					if (gen->parent>=0 && genparts->at(gen->parent)->pdgid==gen->pdgid)
						continue;
					if (processType==kZ) {
						gt->trueGenBosonPt = gen->pt;
						gt->genBosonPt = bound(gen->pt,genBosonPtMin,genBosonPtMax);
						gt->sf_qcdV = hZNLO->Eval(gt->genBosonPt);
						gt->sf_ewkV = hZEWK->Eval(gt->genBosonPt);
						found=true;
					} else if (processType==kW) {
						gt->trueGenBosonPt = gen->pt;
						gt->genBosonPt = bound(gen->pt,genBosonPtMin,genBosonPtMax);
						gt->sf_qcdV = hWNLO->Eval(gt->genBosonPt);
						gt->sf_ewkV = hWEWK->Eval(gt->genBosonPt);
						found=true;
					} else if (processType==kA) {
						// take the highest pT
						if (gen->pt > gt->trueGenBosonPt) {
							gt->trueGenBosonPt = gen->pt;
							gt->genBosonPt = bound(gen->pt,genBosonPtMin,genBosonPtMax);
							gt->sf_qcdV = hANLO->Eval(gt->genBosonPt);
							gt->sf_ewkV = hAEWK->Eval(gt->genBosonPt);
						}
					}
				} // target matches
			}

			tr.TriggerEvent("qcd/ewk SFs");
		}

		//lepton SFs
		gt->sf_lep=1; gt->sf_lepReco=1;
		if (!isData) {
			for (unsigned int iL=0; iL!=TMath::Min(gt->nLooseLep,2); ++iL) {
				PObject *lep = looseLeps.at(iL);
				float pt = lep->pt, eta = lep->eta, aeta = TMath::Abs(eta);
				bool isTight = (iL==0 && gt->looseLep1IsTight) || (iL==1 && gt->looseLep2IsTight);
				PMuon *mu = dynamic_cast<PMuon*>(lep);
				if (mu!=NULL) {
					if (gt->npv<=17) {
						if (isTight)
							gt->sf_lep *= hMuTightLoPU->Eval(aeta,pt);
						else
							gt->sf_lep *= hMuLooseLoPU->Eval(aeta,pt);
						gt->sf_lepReco *= hRecoMuLoPU->Eval(aeta,pt);
					} else {
						if (isTight)
							gt->sf_lep *= hMuTightHiPU->Eval(aeta,pt);
						else
							gt->sf_lep *= hMuLooseHiPU->Eval(aeta,pt);
						gt->sf_lepReco *= hRecoMuHiPU->Eval(aeta,pt);
					}
				} else {
					if (isTight)
						gt->sf_lep *= hEleTight->Eval(eta,pt);
					else
						gt->sf_lep *= hEleVeto->Eval(eta,pt);
					gt->sf_lepReco *= hRecoEle->Eval(eta,pt);
				}
			}

			tr.TriggerEvent("lepton SFs");
		}

		//photon SF
		gt->sf_pho=1;
		if (!isData && gt->nLoosePhoton>0) {
			float pt = gt->loosePho1Pt, eta = gt->loosePho1Eta;
			if (gt->loosePho1IsTight)
				gt->sf_pho = hPho->Eval(eta,pt);

			tr.TriggerEvent("photon SFs");
		}

		// scale and PDF weights, if they exist
		if (!isData) {
			unsigned nSyst = geninfo->mcWeights_syst->size();
			auto *w_syst = geninfo->mcWeights_syst;
			gt->scaleUp = 1; gt->scaleDown = 1;
			gt->pdfUp = 1; gt->pdfDown = 1;
			if (nSyst>8) {
				float norm = w_syst->at(0);
				for (unsigned iS=0; iS!=6; ++iS) {
					gt->scale[iS] = w_syst->at(scaleIndices[iS])/norm;
					gt->scaleUp = max(float(gt->scaleUp),float(gt->scale[iS]));
					gt->scaleDown = min(float(gt->scaleDown),float(gt->scale[iS]));
				}
				// now pdf
				if (nSyst>9) {
					float rms = 0;
					for (unsigned iS=9; iS!=nSyst; ++iS) {
						rms += pow(w_syst->at(iS)/norm-1,2);
					}
					rms /= (nSyst-9);
					rms = sqrt(rms);
					gt->pdfUp = 1+rms;
					gt->pdfDown = 1-rms;
				}
			}
			tr.TriggerEvent("qcd uncertainties");
		}

		gt->Fill();

	} // entry loop

	if (DEBUG) { PDebug("PandaAnalyzer::Run","Done with entry loop"); }

} // Run()

