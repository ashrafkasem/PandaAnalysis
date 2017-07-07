#ifndef ANALYZERUTILS_H
#define ANALYZERUTILS_H

// PandaProd Objects
#include "PandaTree/Objects/interface/Event.h"

// PANDACore
#include "PandaCore/Tools/interface/Common.h"
#include "PandaCore/Tools/interface/DataTools.h"
#include "PandaCore/Tools/interface/JERReader.h"

////////////////////////////////////////////////////////////////////////////////////

inline double TTNLOToNNLO(double pt) {
    double a = 0.1102;
    double b = 0.1566;
    double c = -3.685e-4;
    double d = 1.098;

    return TMath::Min(1.25,
                        a*TMath::Exp(-b*pow(pt,2)+1) + c*pt + d);
}

////////////////////////////////////////////////////////////////////////////////////

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


////////////////////////////////////////////////////////////////////////////////////
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
    ~THCorr() {} // does not own histogram!
    double Eval(double x) {
        if (dim!=1) {
          PError("THCorr1::Eval",
              TString::Format("Trying to access a non-1D histogram (%s)!",h->GetName()));
          return -1;
        }
        return getVal(h,bound(x,lo1,hi1));
    }

    double Eval(double x, double y) {
        if (dim!=2) {
          PError("THCorr1::Eval",
             TString::Format("Trying to access a non-2D histogram (%s)!",h->GetName()));
          return -1;
        }
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

////////////////////////////////////////////////////////////////////////////////////

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

inline bool ElectronIP(double eta, double dxy, double dz) {
  double aeta = fabs(eta);
  if (aeta<1.4442) {
    return (dxy < 0.05 && dz < 0.10) ;
  } else {
    return (dxy < 0.10 && dz < 0.20);
  }
}

////////////////////////////////////////////////////////////////////////////////////


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

////////////////////////////////////////////////////////////////////////////////////

#endif
