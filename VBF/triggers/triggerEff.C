#include <Compression.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include "PandaAnalysis/VBF/triggers/triggerEff.h"

void triggerEff(
TString mode="electrons",
TString inputPath="",
TString batchName="",
bool debug=false
) {
  assert(mode=="electrons" || mode=="muons");
  
  // initialize output tree
  TFile *outputFile = TFile::Open(("vbf_batchTree_triggerEff"+(batchName == "" ? "" : "_"+batchName)+".root"), "RECREATE", "", ROOT::CompressionSettings(ROOT::kZLIB,9));
  bool passMetTriggers;
  float leptonPt, leptonEta, leptonPhi, pfMetNoMu, pfMetNoMuPhi, jet1Pt, jet1Eta, jet1Phi, jet1PtUp, jet1PtDown, jet2Pt, jet2Eta, jet2Phi, jet2PtUp, jet2PtDown;
  TTree *effTree = new TTree("effTree", "effTree");
  effTree->Branch("passMetTriggers" , &passMetTriggers);
  effTree->Branch("leptonPt"        , &leptonPt       );
  effTree->Branch("leptonEta"       , &leptonEta      );
  effTree->Branch("leptonPhi"       , &leptonPhi      );
  effTree->Branch("pfMetNoMu"       , &pfMetNoMu      );
  effTree->Branch("pfMetNoMuPhi"    , &pfMetNoMuPhi   );
  effTree->Branch("jet1Pt"          , &jet1Pt         );
  effTree->Branch("jet1Eta"         , &jet1Eta        );
  effTree->Branch("jet1Phi"         , &jet1Phi        );
  effTree->Branch("jet1PtUp"        , &jet1PtUp       );
  effTree->Branch("jet1PtDown"      , &jet1PtDown     );
  effTree->Branch("jet2Pt"          , &jet2Pt         );
  effTree->Branch("jet2Eta"         , &jet2Eta        );
  effTree->Branch("jet2Phi"         , &jet2Phi        );
  effTree->Branch("jet2PtUp"        , &jet2PtUp       );
  effTree->Branch("jet2PtDown"      , &jet2PtDown     );
  
  TChain input("events");
  if(inputPath!="") {
    input.Add(inputPath);
  } else {
    if(mode=="electrons") {
      input.Add("/mnt/hadoop/scratch/yiiyama/pandaf/006/SingleElectron+Run2017A-PromptReco-v2+MINIAOD/*.root");
      input.Add("/mnt/hadoop/scratch/yiiyama/pandaf/006/SingleElectron+Run2017A-PromptReco-v3+MINIAOD/*.root");
      input.Add("/mnt/hadoop/scratch/yiiyama/pandaf/006/SingleElectron+Run2017B-PromptReco-v1+MINIAOD/*.root");
    } else if(mode=="muons") {
  
    }
  }
  panda::Event event;
  event.setStatus(input, {"!*"});
  if     (mode=="electrons")  event.setAddress(input, {"electrons", "pfMet", "triggers", "chsAK4Jets", "runNumber"});
  else if(mode=="muons"    )  event.setAddress(input, {"muons", "pfMet", "triggers", "chsAK4Jets", "runNumber"});
  long iEntry = 0;
  vector<unsigned> lepTriggerTokens, metTriggerTokens;
  if (mode=="electrons") {
    lepTriggerTokens.push_back( event.registerTrigger("HLT_Ele35_WPTight_Gsf") );
    metTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight") );
    metTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight") );
  } else if(mode=="muons") {
    metTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight") );
    metTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight") );
  }
  
  while (event.getEntry(input, iEntry++) > 0) {
    //bool myTriggerFired = event.triggerFired(triggerToken);
    //if(myTriggerFired) printf("cool man\n");
    bool passedLepTriggers=false;
    for(unsigned i=0; i<lepTriggerTokens.size() && !passedLepTriggers; i++) if(event.triggerFired(lepTriggerTokens[i])) passedLepTriggers=true;
    if(!passedLepTriggers) {  if(debug) printf("failed lepton triggers\n"); continue; }
    
    std::vector<panda::Lepton*> looseLeps, tightLeps; // Fakeable object and tight leptons
    std::vector<panda::Particle*> matchLeps;

    // Find veto and tight electrons
    if(mode=="electrons") for (auto& ele : event.electrons) {
     float pt = ele.pt(); float eta = ele.eta(); float aeta = fabs(eta);
      if (pt<=10 || aeta>2.5 || (aeta>1.4442 && aeta<1.566)) // electron acceptance cuts
        continue;
      if (!ElectronIP(ele.eta(),ele.dxy,ele.dz)) continue;
      if (!ele.veto) continue;
      looseLeps.push_back(&ele); // passes MFO definition
      if(ele.tight) {
        tightLeps.push_back(&ele);
        matchLeps.push_back(&ele);
      }
    }
    
    // Find loose and tight muons
    TVector2 vMETNoMu; vMETNoMu.SetMagPhi(event.pfMet.pt, event.pfMet.phi); // initialize Met no Mu vector as just the PF MET
    if(mode=="muons") for (auto& mu : event.muons) {
      float pt = mu.pt(); float eta = mu.eta(); float aeta = fabs(eta);
      if (pt<=10 || aeta>2.4) continue; //muon acceptance cuts
      if (!mu.loose)  continue;
      if (!MuonIsolation(pt,eta,mu.combIso(),panda::kLoose)) continue;
      looseLeps.push_back(&mu); //passes MFO definition
      if(mu.tight && MuonIsolation(mu.pt(),mu.eta(),mu.combIso(),panda::kTight)) {
        tightLeps.push_back(&mu);
        matchLeps.push_back(&mu);
      }

      // add the muon pT vector to the Met no Mu
      TVector2 vMu; vMu.SetMagPhi(pt,mu.phi());
      vMETNoMu += vMu;
    }
    // only consider events with exactly one lepton which is tight
    if(tightLeps.size()!=1 || looseLeps.size()!=tightLeps.size()) continue;

    // Look for jets 
    panda::JetCollection* jets(0);
    jets = &event.chsAK4Jets;
    vector<panda::Jet*> centralJets;
    TLorentzVector vJet;
    panda::Jet *jet1=0, *jet2=0, *jetUp1=0, *jetUp2=0, *jetDown1=0, *jetDown2=0; 
    for (auto& jet : *jets) {
      // jet cleaning
      if (IsMatched(&matchLeps,0.16,jet.eta(),jet.phi()))
        continue;
      if (!jet.loose)
        continue;

      if (jet.pt()>30) { // consider only central jets here
        if (fabs(jet.eta())<2.4) {
          centralJets.push_back(&jet);
          if (centralJets.size()==1) {
            jet1 = &jet;
            jet1Pt = jet.pt();
            jet1Eta = jet.eta();
            jet1Phi = jet.phi();
          } else if (centralJets.size()==2) {
            jet2 = &jet;
            jet2Pt = jet.pt();
            jet2Eta = jet.eta();
            jet2Phi = jet.phi();
          }
        }
      }

      // do jes variation OUTSIDE of pt>30 check
      if (jet.ptCorrUp>30) {
        if (jet.ptCorrUp > jet1PtUp) {
          if (jetUp1) {
            jetUp2 = jetUp1;
            jet2PtUp  = jet1PtUp;
          }
          jetUp1 = &jet;
          jet1PtUp = jet.ptCorrUp;
        } else if (jet.ptCorrUp > jet2PtUp) {
          jetUp2 = &jet;
          jet2PtUp = jet.ptCorrUp;
        }
      }
      if (jet.ptCorrDown>30) {
        if (jet.ptCorrDown > jet1PtDown) {
          if (jetDown1) {
            jetDown2 = jetDown1;
            jet2PtDown  = jet1PtDown;
          }
          jetDown1 = &jet;
          jet1PtDown = jet.ptCorrDown;
        } else if (jet.ptCorrDown > jet2PtDown) {
          jetDown2 = &jet;
          jet2PtDown = jet.ptCorrDown;
        }
      }
    } //end loop over jets
    
    // Record if the event passed the MET triggers
    passMetTriggers=false;
    for(unsigned i=0; i<metTriggerTokens.size() && !passMetTriggers; i++) if(event.triggerFired(metTriggerTokens[i])) passMetTriggers=true;
    
    leptonPt=tightLeps[0]->pt();
    leptonEta=tightLeps[0]->eta();
    leptonPhi=tightLeps[0]->phi();
    pfMetNoMu=vMETNoMu.Mod();
    pfMetNoMuPhi=vMETNoMu.Phi();
    effTree->Fill();  

    
  }
  outputFile->cd();
  effTree->Write();
  outputFile->Close();
}
