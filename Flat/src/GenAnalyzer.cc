#include "../interface/GenAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define NPARTICLEMAX 64

#define DEBUG 0
using namespace std;

GenAnalyzer::GenAnalyzer() {
  kt = new KFactorTree();
  flags["plots"] = false;
  for (unsigned iH=0; iH!=nPlots; ++iH) {
    hists[iH] = NULL;
  }
}

GenAnalyzer::~GenAnalyzer() {
  for (unsigned iH=0; iH!=nPlots; ++iH) {
    delete hists[iH];
  }
}

void GenAnalyzer::ResetBranches() {
  kt->Reset();
}

void GenAnalyzer::SetOutputPath(TString fOutName, TString fOutFileName) {
  outpath = fOutName;

  fOut = new TFile(fOutName+"/"+fOutFileName,"RECREATE");
  tOut = new TTree("events","events");

  kt->Reset();
  kt->WriteTree(tOut);

}

void GenAnalyzer::AddInput(TTree *t, double xsec, TString label)
{
  if (!t) return;
  InputTree thisTree;
  thisTree.t = t;
  thisTree.xsec = xsec;
  thisTree.label = label;
  tIns.push_back(thisTree);
}

void GenAnalyzer::Terminate() {
  fOut->WriteTObject(tOut);
  fOut->Close();
}

// run
void GenAnalyzer::Run() {

  if (order==kLO) {
    RunLO();
  } else {
    RunNLO();
  }
} // Run()

void GenAnalyzer::RunLO() {

  double weight=0;
  unsigned int nLeps=0, nJets=0;
  int vId;
  float vPt, vEta, vM; 
  int lId[NPARTICLEMAX];
  float lPt[NPARTICLEMAX], lEta[NPARTICLEMAX], lPhi[NPARTICLEMAX], lM[NPARTICLEMAX];
  float jPt[NPARTICLEMAX], jEta[NPARTICLEMAX], jPhi[NPARTICLEMAX], jM[NPARTICLEMAX];

  for (InputTree intree : tIns) {
    TTree *t = intree.t;

    TH1F *hsumw = new TH1F("hsumw","hsumw",1,0,2);
    t->Draw("1>>hsumw","weight");
    double sumw = hsumw->Integral();
    sumw = 1; // I think this is already normalized?
    hsumw->Delete();

    nEvents = t->GetEntries();
    nZero = 0;
    if (lastEvent>=0 && lastEvent<(int)nEvents)
      nEvents = lastEvent;
    if (firstEvent>=0)
      nZero = firstEvent;

    
    // set branch addresses
    t->SetBranchStatus("*",0);
    t->SetBranchStatus("weight",1);
    t->SetBranchStatus("v*",1);
    t->SetBranchStatus("nLeps",1);
    t->SetBranchStatus("l*",1);
    t->SetBranchStatus("nJets",1);
    t->SetBranchStatus("j*",1);

    t->SetBranchAddress("weight",&weight);
    t->SetBranchAddress("nLeps",&nLeps);
    t->SetBranchAddress("nJets",&nJets);
    t->SetBranchAddress("vId",&vId);
    t->SetBranchAddress("vPt",&vPt);
    t->SetBranchAddress("vEta",&vEta);
    t->SetBranchAddress("vM",&vM);
    t->SetBranchAddress("lId",&lId);
    t->SetBranchAddress("lPt",&lPt);
    t->SetBranchAddress("lEta",&lEta);
    t->SetBranchAddress("lPhi",&lPhi);
    // t->SetBranchAddress("lM",&lM);
    t->SetBranchAddress("jPt",&jPt);
    t->SetBranchAddress("jEta",&jEta);
    t->SetBranchAddress("jPhi",&jPhi);
    t->SetBranchAddress("jM",&jM);


    // set up reporters
    unsigned int iE=nZero;
    ProgressReporter pr("GenAnalyzer::Run",&iE,&nEvents,10);
    TimeReporter tr("GenAnalyzer::Run",DEBUG);

    for (iE=nZero; iE!=nEvents; ++iE) {
      pr.Report();
      ResetBranches();

      t->GetEntry(iE);
      kt->vpt = vPt;

      if (kt->vpt<100)
        continue;

      kt->eventNumber =  iE;
      kt->weight = weight*intree.xsec/sumw;
      kt->vid = vId;
      kt->veta = vEta;
      kt->vm = vM;

      TLorentzVector vMET(0,0,0,0); 
      for (unsigned iL=0; iL!=std::min((unsigned)2,nLeps); ++iL) {
        if (iL==0) {
          kt->lep1id = lId[iL];
          kt->lep1pt = lPt[iL];
          kt->lep1eta = lEta[iL];
        } else if (iL==1) {
          kt->lep2id = lId[iL];
          kt->lep2pt = lPt[iL];
          kt->lep2eta = lEta[iL];
        }
        if (lId[iL]%2==0 /*|| TMath::Abs(lEta[iL])>2.4 */) {
          TLorentzVector vNu; vNu.SetPtEtaPhiM(lPt[iL],lEta[iL],lPhi[iL],0);
          vMET += vNu;
        }
      }

      kt->met = vMET.Pt();

      TLorentzVector vj1, vj2;
      kt->njet = 0;
      for (unsigned iJ=0; iJ!=nJets; ++iJ) {
        kt->njet++;
        kt->ht += jPt[iJ];
        if (kt->njet==1) {
          kt->jet1pt = jPt[iJ];
          kt->jet1eta = jEta[iJ];
          vj1.SetPtEtaPhiM(jPt[iJ],jEta[iJ],jPhi[iJ],jM[iJ]);
        } else if (kt->njet==2) {
          kt->jet2pt = jPt[iJ];
          kt->jet2eta = jEta[iJ];
          vj2.SetPtEtaPhiM(jPt[iJ],jEta[iJ],jPhi[iJ],jM[iJ]);
        }
      }

      kt->jjdeta = TMath::Abs(vj1.Eta()-vj2.Eta());
      kt->jjdphi = vj1.DeltaPhi(vj2);
      kt->mjj    = (vj1+vj2).M();

      tOut->Fill();
    }
  }
}

void GenAnalyzer::RunNLO() {

  struct Particles {
    unsigned size;
    float pt[NPARTICLEMAX];
    float eta[NPARTICLEMAX];
    float phi[NPARTICLEMAX];
    float mass[NPARTICLEMAX];
  };

  double weight;
  Particles electrons;
  Particles muons;
  Particles neutrinos;
  Particles photons;
  Particles jets;
  float met_pt, met_phi;

  for (InputTree intree : tIns) {
    TTree *t = intree.t;

    TH1F *hsumw = new TH1F("hsumw","hsumw",1,0,2);
    t->Draw("1>>hsumw","weight");
    double sumw = hsumw->Integral();
    hsumw->Delete();

    nEvents = t->GetEntries();
    nZero = 0;
    if (lastEvent>=0 && lastEvent<(int)nEvents)
      nEvents = lastEvent;
    if (firstEvent>=0)
      nZero = firstEvent;
    
    // set branch addresses
    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("weight", 1);
    t->SetBranchStatus("electrons.*", 1);
    t->SetBranchStatus("muons.*", 1);
    t->SetBranchStatus("neutrinos.*", 1);
    t->SetBranchStatus("jets.*", 1);
    t->SetBranchStatus("met.*", 1);

    t->SetBranchAddress("weight", &weight);
    t->SetBranchAddress("electrons.size", &electrons.size);
    t->SetBranchAddress("electrons.pt", electrons.pt);
    t->SetBranchAddress("electrons.eta", electrons.eta);
    t->SetBranchAddress("electrons.phi", electrons.phi);
    t->SetBranchAddress("muons.size", &muons.size);
    t->SetBranchAddress("muons.pt", muons.pt);
    t->SetBranchAddress("muons.eta", muons.eta);
    t->SetBranchAddress("muons.phi", muons.phi);
    t->SetBranchAddress("neutrinos.size", &neutrinos.size);
    t->SetBranchAddress("neutrinos.pt", neutrinos.pt);
    t->SetBranchAddress("neutrinos.eta", neutrinos.eta);
    t->SetBranchAddress("neutrinos.phi", neutrinos.phi);
    t->SetBranchAddress("photons.size", &photons.size);
    t->SetBranchAddress("photons.pt", photons.pt);
    t->SetBranchAddress("photons.eta", photons.eta);
    t->SetBranchAddress("photons.phi", photons.phi);
    t->SetBranchAddress("jets.size", &jets.size);
    t->SetBranchAddress("jets.pt", jets.pt);
    t->SetBranchAddress("jets.eta", jets.eta);
    t->SetBranchAddress("jets.phi", jets.phi);
    t->SetBranchAddress("jets.mass", jets.mass);
    t->SetBranchAddress("met.pt", &met_pt);
    t->SetBranchAddress("met.phi", &met_phi);

    // set up reporters
    unsigned int iE=nZero;
    ProgressReporter pr("GenAnalyzer::Run",&iE,&nEvents,10);
    TimeReporter tr("GenAnalyzer::Run",DEBUG);

    for (iE=nZero; iE!=nEvents; ++iE) {
      pr.Report();
      ResetBranches();

      t->GetEntry(iE);

      kt->eventNumber =  iE;
      kt->weight = weight*intree.xsec/sumw;

      kt->met = met_pt;

      TLorentzVector vV(0,0,0,0);

      if (processType==kZ) {
        kt->vid = 23;
        // neutrinos only for NLO Zvv
        TLorentzVector vNu;
        for (unsigned iNu=0; iNu!=neutrinos.size; ++iNu) {
          vNu.SetPtEtaPhiM(neutrinos.pt[iNu],neutrinos.eta[iNu],neutrinos.phi[iNu],0);          
          vV += vNu;
        }
        //kt->vpt = vV.Pt();
        //if (kt->vpt<100)
        //  continue;
        kt->vpt = kt->met;
        if (kt->met<100) // sometimes neutrinos are missing ???
          continue;
      } else {
        if (electrons.size==0 && muons.size==0)
          continue;

        TLorentzVector vL, vNu, vNu_tmp;
        // first search for an e-nu_e pair
        for (unsigned iE=0; iE!=electrons.size; ++iE) {
          bool decayLike = false;
          for (unsigned iNu=0; iNu!=neutrinos.size; ++iNu) {
            if (DeltaR2(electrons.eta[iE],electrons.phi[iE],neutrinos.eta[iNu],neutrinos.phi[iNu])<0.04) {
              decayLike = true;
              break;
            }
            vNu_tmp.SetPtEtaPhiM(neutrinos.pt[iNu],neutrinos.eta[iNu],neutrinos.phi[iNu],0);
            break;
          }
          if (decayLike)
            continue;
          if (electrons.pt[iE] > vL.Pt()) {
            vL.SetPtEtaPhiM(electrons.pt[iE],electrons.eta[iE],electrons.phi[iE],0.000511);
            vNu = vNu_tmp;
            kt->lep1id = 11; // do we care about the sign?
            kt->lep1pt = electrons.pt[iE];
            kt->lep1eta = electrons.eta[iE];
          }
        }
        // first search for an mu-nu_mu pair
        for (unsigned iM=0; iM!=muons.size; ++iM) {
          bool decayLike = false;
          for (unsigned iNu=0; iNu!=neutrinos.size; ++iNu) {
            if (DeltaR2(muons.eta[iM],muons.phi[iM],neutrinos.eta[iNu],neutrinos.phi[iNu])<0.04) {
              decayLike = true;
              break;
            }
            vNu_tmp.SetPtEtaPhiM(neutrinos.pt[iNu],neutrinos.eta[iNu],neutrinos.phi[iNu],0);
            break;
          }
          if (decayLike)
            continue;
          if (muons.pt[iM] > vL.Pt()) {
            vL.SetPtEtaPhiM(muons.pt[iM],muons.eta[iM],muons.phi[iM],0.105);
            vNu = vNu_tmp;
            kt->lep1id = 13; // do we care about the sign?
            kt->lep1pt = muons.pt[iM];
            kt->lep1eta = muons.eta[iM];
          }
        }
        vV = vL+vNu;
        kt->vpt = vV.Pt();
        if (kt->vpt<100)
          continue;
      }

      kt->veta = vV.Eta();
      kt->vm = vV.M();

      TLorentzVector vj1, vj2;
      kt->njet = 0;
      for (unsigned iJ=0; iJ!=jets.size; ++iJ) {
        kt->njet++;
        kt->ht += jets.pt[iJ];
        if (kt->njet==1) {
          kt->jet1pt = jets.pt[iJ];
          kt->jet1eta = jets.eta[iJ];
          vj1.SetPtEtaPhiM(jets.pt[iJ],jets.eta[iJ],jets.phi[iJ],jets.mass[iJ]);
        } else if (kt->njet==2) {
          kt->jet2pt = jets.pt[iJ];
          kt->jet2eta = jets.eta[iJ];
          vj2.SetPtEtaPhiM(jets.pt[iJ],jets.eta[iJ],jets.phi[iJ],jets.mass[iJ]);
        }
      }

      kt->jjdeta = TMath::Abs(vj1.Eta()-vj2.Eta());
      kt->jjdphi = vj1.DeltaPhi(vj2);
      kt->mjj    = (vj1+vj2).M();

      tOut->Fill();
    }
  }
}
