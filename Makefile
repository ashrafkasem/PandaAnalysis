# only generates ROOT dictionaries! compilation is left to SCRAM
PACKAGE=$(CMSSW_BASE)/src/PandaAnalysis/
INC=-I/cvmfs/cms.cern.ch/$(SCRAM_ARCH)/external/fastjet-contrib/1.020/include/ -I/cvmfs/cms.cern.ch/$(SCRAM_ARCH)/external/fastjet/3.1.0/include/

.PHONY: clean all

all: Flat/src/dictFlat.cc 

clean:
	rm -f */src/dict*.cc */src/dict*pcm

Flat/src/dictFlat.cc: $(wildcard $(PACKAGE)/Flat/interface/*.h) $(PACKAGE)/Flat/LinkDef.h
	rootcling -f $(PACKAGE)/Flat/src/dictFlat.cc $(CMSSW_BASE)/src/PandaAnalysis/Flat/interface/*.h $(INC) $(CMSSW_BASE)/src/PandaAnalysis/Flat/LinkDef.h
	mkdir -p $(CMSSW_BASE)/lib/$(SCRAM_ARCH)/
	cp Flat/src/dictFlat_rdict.pcm $(CMSSW_BASE)/lib/$(SCRAM_ARCH)/


# DO NOT DELETE
