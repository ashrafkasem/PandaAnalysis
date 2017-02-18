# only generates ROOT dictionaries! compilation is left to SCRAM
PACKAGE=$(CMSSW_BASE)/src/PandaAnalysis/

.PHONY: clean all

all: Flat/src/dictFlat.cc 

clean:
	rm -f */src/dict*.cc */src/dict*pcm

Flat/src/dictFlat.cc: $(wildcard $(PACKAGE)/Flat/interface/*.h) $(PACKAGE)/Flat/LinkDef.h
	rootcling -f $(PACKAGE)/Flat/src/dictFlat.cc $(CMSSW_BASE)/src/PandaAnalysis/Flat/interface/*.h $(CMSSW_BASE)/src/PandaAnalysis/Flat/LinkDef.h
	mkdir -p $(CMSSW_BASE)/lib/$(SCRAM_ARCH)/
	cp Flat/src/dictFlat_rdict.pcm $(CMSSW_BASE)/lib/$(SCRAM_ARCH)/
