#!/usr/bin/env python

from re import sub
from sys import argv,exit
from os import system,getenv,path
from time import clock
import subprocess

which = int(argv[1])
sname = argv[0]
if len(argv)>2:
	nper = int(argv[2])
else:
	nper = 1
argv=[]

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import *
import PandaAnalysis.Tagging.cfg_v8 as tagcfg

if __name__ == "__main__":
	
	Load('PandaAnalysisFlat','PandaAnalyzer')

	def fn(shortName,longName,isData,xsec):
		# first we do some I/O stuff
		outdir = 'XXXX'
		outfilename = shortName+'.root'
		PInfo(sname,'Output: %s/%s'%(outdir,outfilename))
		if path.isfile('%s/%s'%(outdir,outfilename)):
			PWarning(sname,"Found %s, skipping!"%outfilename)
			return;

		start=clock()

		eosPath = 'root://eoscms.cern.ch//store/user/snarayan'
		eosEXOPath = 'root://eoscms.cern.ch//store/group/phys_exotica'
		cernboxPath = 'root://eosuser//eos/user/s/snarayan'
		cernboxBPath = 'root://eosuser//eos/user/b/bmaier'
		fullPath = sub(r'\${CERNBOXB}',cernboxBPath,
				sub(r'\${CERNBOX}',cernboxPath,
					sub(r'\${EOS}',eosPath,
						sub(r'\${EOS2}',eosEXOPath,longName))))
		PInfo(sname,fullPath)

		# xrdcp if remote, copy if local
		if 'root://' in fullPath:
			system('xrdcp %s input.root'%fullPath)
		else:
			system('cp %s input.root'%fullPath)

		# now we instantiate and configure the analyzer
		skimmer = root.PandaAnalyzer()
		skimmer.isData=isData
		skimmer.SetFlag('firstGen',False)
		# skimmer.SetFlag('puppi',False)
		skimmer.SetPreselectionBit(root.PandaAnalyzer.kRecoil)
		#skimmer.SetPreselectionBit(root.PandaAnalyzer.kMonotop)
		processType=root.PandaAnalyzer.kNone
		if not isData:
			if 'ST_' in fullPath:
				processType=root.PandaAnalyzer.kTop
			elif 'ZJets' in fullPath or 'DY' in fullPath:
				processType=root.PandaAnalyzer.kZ
			elif 'WJets' in fullPath:
				processType=root.PandaAnalyzer.kW
			elif 'GJets' in fullPath:
				processType=root.PandaAnalyzer.kA
			elif 'TTJets' in fullPath or 'TT_' in fullPath:
				processType=root.PandaAnalyzer.kTT
		skimmer.processType=processType 

		# read the inputs
		try:
			fin = root.TFile.Open('input.root')
			tree = fin.FindObjectAny("events")
			infotree = fin.FindObjectAny("all")
		except:
			exit(2) # file open error => xrootd?

		skimmer.SetDataDir(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/')
		skimmer.SetOutputFile('output.root')
		skimmer.Init(tree,infotree)

		# run and save output
		skimmer.Run()
		skimmer.Terminate()

		# now run the BDT
		Load('Learning','TMVABranchAdder')
		ba = root.TMVABranchAdder()
		ba.treename = 'events'
		ba.defaultValue = -1.2
		ba.presel = 'fj1ECFN_2_4_20>0'
		for v in tagcfg.variables:
			ba.AddVariable(v[0],v[2])
		for v in tagcfg.formulae:
			ba.AddFormula(v[0],v[2])
		for s in tagcfg.spectators:
			ba.AddSpectator(s[0])
		ba.BookMVA('top_ecf_bdt',getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/trainings/top_ecfbdt_v8_BDT.weights.xml')
		ba.RunFile('output.root')

		# stageout
		mvargs = 'mv $PWD/output.root %s/%s'%(outdir,outfilename)
		PInfo(sname,mvargs)
		system(mvargs)
		system('rm input.root')

		PInfo(sname,'finished in %f'%(clock()-start)); start=clock()

	cfg = open('local.cfg')
	lines = list(cfg)
	lines_ = lines[which*nper:min(len(lines),(which+1)*nper)]
	PDebug(sname,'%i ->%i'%(which*nper,min(len(lines),(which+1)*nper)))
	PDebug(sname,'%i %i'%(len(lines),len(lines_)))
	del lines
	cfg.close()

	for line in lines_:
		ll = line.split()
		shortname = ll[0]
		isData = (ll[1]!="MC")
		xsec = float(ll[2])
		longname = ll[3]
		fn(shortname,longname,isData,xsec)

	exit(0)
