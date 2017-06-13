#ifndef GENERALLEPTONICTREE_H
#define GENERALLEPTONICTREE_H 

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TString.h"
#include "genericTree.h"
#include <map>

class GeneralLeptonicTree : public genericTree {
    public:
      // public objects

      enum BTagShift {
        bCent=0,
        bBUp,
        bBDown,
        bMUp,
        bMDown,
        bNShift
      };
      enum BTagJet {
        bJet=0,
        bSubJet,
        bNJet
      };
      enum BTagTags {
        b0=0,
        b1,
        b2,
        bGT0,
        bNTags
      };
        
      struct BTagParams {
        BTagJet jet;
        BTagTags tag;
        BTagShift shift=(BTagShift)0;
        bool operator==(const BTagParams &o) const {
          return jet==o.jet && tag==o.tag && shift==o.shift;
        }
        bool operator<(const BTagParams &o) const {
          return ( jet<o.jet ||
                   (jet==o.jet && tag<o.tag) ||
                   (jet==o.jet && tag==o.tag && shift<o.shift) ); 
        }
        bool operator>(const BTagParams &o) const {
          return ! operator<(o);
        }
      };
        
    private:
        std::vector<BTagParams> btagParams;

        TString makeBTagSFString(BTagParams p) {
          TString s = "sf_";
          if (p.jet==bSubJet)
            s += "sj";
          s += "btag";
          switch (p.tag) {
            case b0:
              s += "0"; break;
            case b1:
              s += "1"; break;
            case b2:
              s += "2"; break;
            case bGT0:
              s += "GT0"; break;
            default:
              break;
          }
          if (p.shift==bCent)
            return s;
          switch (p.shift) {
            case bBUp:
              s += "BUp"; break;
            case bBDown:
              s += "BDown"; break;
            case bMUp:
              s += "MUp"; break;
            case bMDown:
              s += "MDown"; break;
            default: break;
          }
          return s;
        }
    public:
      GeneralLeptonicTree();
      ~GeneralLeptonicTree();
      void WriteTree(TTree *t);
      void Fill() { treePtr->Fill(); }
      void SetBranchStatus(const char *bname, bool status, UInt_t *ret=0) 
      { 
        treePtr->SetBranchStatus(bname,status,ret); 
      }
      void Reset();

      // public config

//STARTCUSTOMDEF
      std::map<BTagParams,float> sf_btags;
      std::map<TString,float> signal_weights;

      float scale[6];
//ENDCUSTOMDEF
    int runNumber = -1;
    int lumiNumber = -1;
    ULong64_t eventNumber = -1;
    int npv = -1;
    int pu = -1;
    float mcWeight = 1;
    int trigger = 0;
    int metFilter = 0;
    int egmFilter = 0;

    int nLooseLep = 0;
    int looseGenLep1PdgId = 0;
    int looseGenLep2PdgId = 0;
    int looseGenLep3PdgId = 0;
    int looseGenLep4PdgId = 0;
    int looseLep1PdgId = -1;
    int looseLep2PdgId = -1;
    int looseLep3PdgId = -1;
    int looseLep4PdgId = -1;
    int looseLep1SelBit = 0;
    int looseLep2SelBit = 0;
    int looseLep3SelBit = 0;
    int looseLep4SelBit = 0;
    float looseLep1Pt = -1;
    float looseLep2Pt = -1;
    float looseLep3Pt = -1;
    float looseLep4Pt = -1;
    float looseLep1Eta = -1;
    float looseLep2Eta = -1;
    float looseLep3Eta = -1;
    float looseLep4Eta = -1;
    float looseLep1Phi = -1;
    float looseLep2Phi = -1;
    float looseLep3Phi = -1;
    float looseLep4Phi = -1;

    int nJet = 0;
    int jetNLBtags = 0;
    int jetNMBtags = 0;
    int jetNTBtags = 0;

    float jet1Pt    = -1;
    float jet2Pt    = -1;
    float jet3Pt    = -1;
    float jet4Pt    = -1;
    float jet1Eta   = -1;
    float jet2Eta   = -1;
    float jet3Eta   = -1;
    float jet4Eta   = -1;
    float jet1Phi   = -1;
    float jet2Phi   = -1;
    float jet3Phi   = -1;
    float jet4Phi   = -1;
    float jet1BTag  = -1;
    float jet2BTag  = -1;
    float jet3BTag  = -1;
    float jet4BTag  = -1;
    float jet1GenPt = -1;
    float jet2GenPt = -1;
    float jet3GenPt = -1;
    float jet4GenPt = -1;
    int jet1Flav    = -1;
    int jet2Flav    = -1;
    int jet3Flav    = -1;
    int jet4Flav    = -1;
    int jet1SelBit  = 0;
    int jet2SelBit  = 0;
    int jet3SelBit  = 0;
    int jet4SelBit  = 0;

    float jet1PtUp    = -1;
    float jet2PtUp    = -1;
    float jet3PtUp    = -1;
    float jet4PtUp    = -1;
    float jet1PtDown  = -1;
    float jet2PtDown  = -1;
    float jet3PtDown  = -1;
    float jet4PtDown  = -1;
    float jet1EtaUp   = -1;
    float jet2EtaUp   = -1;
    float jet3EtaUp   = -1;
    float jet4EtaUp   = -1;
    float jet1EtaDown = -1;
    float jet2EtaDown = -1;
    float jet3EtaDown = -1;
    float jet4EtaDown = -1;

    float pfmet = -1;
    float pfmetphi = -1;
    float pfmetRaw = -1;
    float pfmetUp = -1;
    float pfmetDown = -1;
    float pfmetnomu = -1;
    float puppimet = -1;
    float puppimetphi = -1;
    float calomet = -1;
    float calometphi = -1;
    float trkmet = -1;
    float trkmetphi = -1;
    float dphipfmet = -1;
    float dphipuppimet = -1;

    float genLep1Pt = -1;
    float genLep2Pt = -1;
    float genLep1Eta = -1;
    float genLep2Eta = -1;
    float genLep1Phi = -1;
    float genLep2Phi = -1;
    int genLep1PdgId = -1;
    int genLep2PdgId = -1;

    int nTau = 0;
    float pdfUp = 1;
    float pdfDown = 1;

    int nLoosePhoton = 0;
    float loosePho1Pt  = -1;
    float loosePho1Eta = -1;
    float loosePho1Phi = -1;
    
    float sf_pu = 1.0;
    float sf_puUp = 1.0;
    float sf_puDown = 1.0;
    float sf_zz = 1.0;
    float sf_zzUnc = 1.0;
    float sf_wz = 1.0;
    float sf_zh = 1.0;
    float sf_zhUp = 1.0;
    float sf_zhDown = 1.0;
    float sf_tt = 1.0;

    float sf_trk1 = 1.0;
    float sf_trk2 = 1.0;
    float sf_trk3 = 1.0;
    float sf_trk4 = 1.0;
    float sf_loose1 = 1.0;
    float sf_loose2 = 1.0;
    float sf_loose3 = 1.0;
    float sf_loose4 = 1.0;
    float sf_medium1 = 1.0;
    float sf_medium2 = 1.0;
    float sf_medium3 = 1.0;
    float sf_medium4 = 1.0;
    float sf_tight1 = 1.0;
    float sf_tight2 = 1.0;
    float sf_tight3 = 1.0;
    float sf_tight4 = 1.0;

};

#endif
