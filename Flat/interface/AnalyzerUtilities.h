#ifndef ANALYZERUTILS_H
#define ANALYZERUTILS_H

// PandaProd Objects
#include "PandaTree/Objects/interface/Event.h"

// PANDACore
#include "PandaCore/Tools/interface/Common.h"
#include "PandaCore/Tools/interface/DataTools.h"
#include "PandaCore/Tools/interface/JERReader.h"

// fastjet
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/MeasureDefinition.hh"

////////////////////////////////////////////////////////////////////////////////////
typedef std::vector<fastjet::PseudoJet> VPseudoJet;

inline VPseudoJet ConvertPFCands(std::vector<const panda::PFCand*> &incoll, bool puppi, double minPt=0.001) {
  VPseudoJet vpj;
  vpj.reserve(incoll.size());
  for (auto *incand : incoll) {
    double factor = puppi ? incand->puppiW() : 1;
    if (factor*incand->pt()<minPt)
      continue;
    vpj.emplace_back(factor*incand->px(),factor*incand->py(),
                     factor*incand->pz(),factor*incand->e());
  }
  return vpj;
}

inline VPseudoJet ConvertPFCands(panda::RefVector<panda::PFCand> &incoll, bool puppi, double minPt=0.001) {
  std::vector<const panda::PFCand*> outcoll;
  outcoll.reserve(incoll.size());
  for (auto incand : incoll)
    outcoll.push_back(incand.get());

  return ConvertPFCands(outcoll, puppi, minPt);
}

inline VPseudoJet ConvertPFCands(panda::PFCandCollection &incoll, bool puppi, double minPt=0.001) {
  std::vector<const panda::PFCand*> outcoll;
  outcoll.reserve(incoll.size());
  for (auto &incand : incoll)
    outcoll.push_back(&incand);

  return ConvertPFCands(outcoll, puppi, minPt);
}

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

class Analysis {
public:
  Analysis(TString name_ = "") { name = name_; }
  ~Analysis() {}
  TString name;
  bool rerunJES = false;
  bool varyJES = false;
  bool complicatedLeptons = false;
  bool vbf = false;
  bool recoil = true;
  bool fatjet = true;
  bool monoh = false;
  bool recluster = false;
  bool genOnly = false;
  bool btagSFs = true;
  bool firstGen = true;
  bool puppi_jets = true;
  bool ak8 = false;
  bool reclusterGen = false;
};

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
class TriggerHandler {  
public:
  TriggerHandler() {};
  ~TriggerHandler() {};
  void addTriggers(std::vector<TString> paths) { 
    for (auto &path : paths) {
      paths.push_back(path); 
      indices.push_back(-1); 
    }
  }
  void registerTrigger(unsigned my_idx, int panda_idx) { indices[my_idx] = panda_idx; }
  std::vector<int> indices;
  std::vector<TString> paths;
};


////////////////////////////////////////////////////////////////////////////////////
template <typename T>
class TCorr {
public:
  TCorr(T*) {} 
  virtual ~TCorr() {}
  virtual double Eval(double x)=0; 
protected:
  T *h=0;
};


class TF1Corr : public TCorr<TF1> {
public:
  TF1Corr(TF1 *f_):
    TCorr(f_) 
  {
//    f_->SetDirectory(0);
    h = f_;
  }
  ~TF1Corr() {} 
  double Eval(double x) {
    return h->Eval(x);
  }

  TF1 *GetFunc() { return h; }
};


template <typename T>
class THCorr : public TCorr<T> {
public:
  // wrapper around TH* to do corrections
  THCorr(T *h_):
    TCorr<T>(h_)
  {
//    h_->SetDirectory(0);
    this->h = h_;
    dim = this->h->GetDimension();
    TAxis *thurn = this->h->GetXaxis(); 
    lo1 = thurn->GetBinCenter(1);
    hi1 = thurn->GetBinCenter(thurn->GetNbins());
    if (dim>1) {
      TAxis *taxis = this->h->GetYaxis();
      lo2 = taxis->GetBinCenter(1);
      hi2 = taxis->GetBinCenter(taxis->GetNbins());
    }
  }
  ~THCorr() {} // does not own histogram!
  double Eval(double x) {
    if (dim!=1) {
      PError("THCorr1::Eval",
        TString::Format("Trying to access a non-1D histogram (%s)!",this->h->GetName()));
      return -1;
    }
    return getVal(this->h,bound(x,lo1,hi1));
  }

  double Eval(double x, double y) {
    if (dim!=2) {
      PError("THCorr1::Eval",
       TString::Format("Trying to access a non-2D histogram (%s)!",this->h->GetName()));
      return -1;
    }
    return getVal(this->h,bound(x,lo1,hi1),bound(y,lo2,hi2));
  }

  double Error(double x) {
    if (dim!=1) {
      PError("THCorr1::Eval",
        TString::Format("Trying to access a non-1D histogram (%s)!",this->h->GetName()));
      return -1;
    }
    return getError(this->h,bound(x,lo1,hi1));
  }

  double Error(double x, double y) {
    if (dim!=2) {
      PError("THCorr1::Eval",
       TString::Format("Trying to access a non-2D histogram (%s)!",this->h->GetName()));
      return -1;
    }
    return getError(this->h,bound(x,lo1,hi1),bound(y,lo2,hi2));
  }

  T *GetHist() { return this->h; }

private:
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
