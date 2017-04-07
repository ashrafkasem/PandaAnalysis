#include "PandaAnalysis/Flat/interface/PandaAnalyzer.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"

#include <iostream>

int
main(int argc, char const* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: pana input output" << std::endl;
    return 1;
  }

  TFile* source(TFile::Open(argv[1]));
  auto* tree(static_cast<TTree*>(source->Get("events")));

  PandaAnalyzer skimmer;

  skimmer.isData = false;
  skimmer.SetFlag("puppi",true);
  skimmer.SetFlag("fatjet",true);
  skimmer.SetFlag("firstGen",false);
  skimmer.processType = PandaAnalyzer::kW;

  skimmer.SetDataDir(TString(gSystem->Getenv("CMSSW_BASE")) + "/src/PandaAnalysis/data/");
  skimmer.SetOutputFile(argv[2]);
  skimmer.Init(tree, 0);
  skimmer.Run();
  skimmer.Terminate();

  return 0;
}
