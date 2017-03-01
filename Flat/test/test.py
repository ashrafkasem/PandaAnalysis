#!/usr/bin/env python

from re import sub
from sys import argv,exit
from os import system,getenv
import json

torun = argv[1]
argv = []

import ROOT as root
from PandaCore.Tools.Load import *

if __name__ == "__main__":

	Load('PandaAnalysisFlat','PandaAnalyzer')
	
	def fn(fullPath):

		skimmer = root.PandaAnalyzer()
	 
#		skimmer.firstEvent=0
#		skimmer.lastEvent=10
		skimmer.isData=False
		skimmer.SetFlag('puppi',True)
		skimmer.SetFlag('fatjet',True)
		skimmer.SetFlag('firstGen',False)
		skimmer.SetFlag('applyJSON',False)
		skimmer.SetFlag('monohiggs',True)
		if skimmer.isData and False:
			with open(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/certs/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt') as jsonFile:
				payload = json.load(jsonFile)
				for run,lumis in payload.iteritems():
					for l in lumis:
						skimmer.AddGoodLumiRange(int(run),l[0],l[1])
		skimmer.processType = root.PandaAnalyzer.kTT
#		skimmer.SetPreselectionBit(root.PandaAnalyzer.kMonotop)
		fin = root.TFile.Open(fullPath)

		print fullPath
		print fin

		tree = fin.FindObjectAny("events")
		infotree = fin.FindObjectAny("all")
		print tree,infotree

		skimmer.SetDataDir(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/')
		skimmer.SetOutputFile('testskim.root')
		skimmer.Init(tree,infotree)

		skimmer.Run()
		print 'done running'
		skimmer.Terminate()
		print 'done terminating'

	fn(torun) 
