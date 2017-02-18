#ifndef ANALYZERUTILS_H
#define ANALYZERUTILS_H

// PandaProd Objects
#include "PandaTree/Objects/interface/Event.h"

// PANDACore
#include "PandaCore/Tools/interface/Common.h"
#include "PandaCore/Tools/interface/DataTools.h"
#include "PandaCore/Tools/interface/JERReader.h"

inline double TTNLOToNNLO(double pt) {
    double a = 0.1102;
    double b = 0.1566;
    double c = -3.685e-4;
    double d = 1.098;

    return TMath::Min(1.25,
                        a*TMath::Exp(-b*pow(pt,2)+1) + c*pt + d);
}

class LumiRange {
public:
    LumiRange(int l0_,int l1_):
        l0(l0_),
        l1(l1_)
     { }
    ~LumiRange() {}
    bool Contains(int l) {
        return l0<=l && l<=l1;
    }
private:
    int l0, l1;
};


template <typename T>
class THCorr {
public:
    // wrapper around TH* to do corrections
    THCorr(T *h_) {
        h = h_;
        dim = h->GetDimension();
        TAxis *thurn = h->GetXaxis(); 
        lo1 = thurn->GetBinCenter(1);
        hi1 = thurn->GetBinCenter(thurn->GetNbins());
        if (dim>1) {
            TAxis *taxis = h->GetYaxis();
            lo2 = taxis->GetBinCenter(1);
            hi2 = taxis->GetBinCenter(taxis->GetNbins());
        }
    }
    ~THCorr() {}
    double Eval(double x) {
        if (dim!=1)
            return -1;
        return getVal(h,bound(x,lo1,hi1));
    }

    double Eval(double x, double y) {
        if (dim!=2)
            return -1;
        return getVal(h,bound(x,lo1,hi1),bound(y,lo2,hi2));
    }

    T *GetHist() { return h; }    

private:
    T *h;
    int dim;
    double lo1, lo2, hi1, hi2;
};

typedef THCorr<TH1D> THCorr1;
typedef THCorr<TH2D> THCorr2;

namespace panda {
  enum IDWorkingPoint {
    kVeto,
    kLoose,
    kMedium,
    kTight,
    nIDWorkingPoints
  };
}

inline bool MuonIsolation(double pt, double eta, double iso, panda::IDWorkingPoint isoType) {
    float maxIso=0;
    maxIso = (isoType == panda::kTight) ? 0.15 : 0.25;
    return (iso < pt*maxIso);
}

inline bool ElectronIsolation(double pt, double eta, double iso, panda::IDWorkingPoint isoType) {
    float maxIso=0;
    float abseta = TMath::Abs(eta);
    switch (isoType) {
        case panda::kVeto:
            maxIso = (abseta<=1.479) ? 0.126 : 0.144;
            break;
        case panda::kLoose:
            maxIso = (abseta<=1.479) ? 0.0893 : 0.121;
            break;
        case panda::kMedium:
            maxIso = (abseta<=1.479) ? 0.0766 : 0.0678;
            break;
        case panda::kTight:
            maxIso = (abseta<=1.479) ? 0.0354 : 0.0646;
            break;
        default:
            break;
    }
    return (iso < pt*maxIso);
}

inline bool IsMatched(std::vector<panda::Particle*>*objects,
               double deltaR2, double eta, double phi) {
  for (auto *x : *objects) {
    if (x->pt()>0) {
      if ( DeltaR2(x->eta(),x->phi(),eta,phi) < deltaR2 )
        return true;
    }
  }
  return false;
}

class btagcand {
    public:
        btagcand(unsigned int i, int f,double e,double cent,double up,double down) {
            idx = i;
            flav = f;
            eff = e;
            sf = cent;
            sfup = up;
            sfdown = down;
        }
        ~btagcand() { }
        int flav, idx;
        double eff, sf, sfup, sfdown;
};


inline void EvalBtagSF(std::vector<btagcand> cands, std::vector<double> sfs, float &sf0, float &sf1) {
    sf0 = 1; sf1 = 1;
    float prob_mc0=1, prob_data0=1;
    float prob_mc1=0, prob_data1=0;
    unsigned int nC = cands.size();

    for (unsigned int iC=0; iC!=nC; ++iC) {
        double sf_i = sfs[iC];
        double eff_i = cands[iC].eff;
        prob_mc0 *= (1-eff_i);
        prob_data0 *= (1-sf_i*eff_i);
        float tmp_mc1=1, tmp_data1=1;
        for (unsigned int jC=0; jC!=nC; ++jC) {
            if (iC==jC) continue;
            double sf_j = sfs[jC];
            double eff_j = cands[jC].eff;
            tmp_mc1 *= (1-eff_j);
            tmp_data1 *= (1-eff_j*sf_j);
        }
        prob_mc1 += eff_i * tmp_mc1;
        prob_data1 += eff_i * sf_i * tmp_data1;
    }
    
    if (nC>0) {
        sf0 = prob_data0/prob_mc0;
        sf1 = prob_data1/prob_mc1;
    }
}

inline void EvalBtagSF(std::vector<btagcand> cands, std::vector<double> sfs, float &sf0, float &sf1, float &sf2) { 
    EvalBtagSF(cands,sfs,sf0,sf1); // get 0,1

    sf2=1;
    float prob_mc2=0, prob_data2=0;
    unsigned int nC = cands.size();

    for (unsigned int iC=0; iC!=nC; ++iC) {
        double sf_i = sfs[iC], eff_i = cands[iC].eff;
        for (unsigned int jC=iC+1; jC!=nC; ++jC) {
            double sf_j = sfs[jC], eff_j = cands[jC].eff;
            float tmp_mc2=1, tmp_data2=1;
            for (unsigned int kC=0; kC!=nC; ++kC) {
                if (kC==iC || kC==jC) continue;
                double sf_k = sfs[kC], eff_k = cands[kC].eff;
                tmp_mc2 *= (1-eff_k);
                tmp_data2 *= (1-eff_k*sf_k);
            }
            prob_mc2 += eff_i * eff_j * tmp_mc2;
            prob_data2 += eff_i * sf_i * eff_j * sf_j * tmp_data2;
        }
    }

    if (nC>1) {
        sf2 = prob_data2/prob_mc2;
    }
}


inline TString makeECFString(int order, int N, float beta) {
    return TString::Format("ECFN_%i_%i_%.2i",order,N,int(10*beta));
}

#endif
