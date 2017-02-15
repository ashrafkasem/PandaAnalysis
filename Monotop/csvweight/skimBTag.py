#!/usr/bin/env python

from sys import argv,exit
import argparse

parser = argparse.ArgumentParser(description='skip stuff')
parser.add_argument('--infile',type=str)
parser.add_argument('--outfile',type=str)
args = parser.parse_args()
sname = argv[0]
argv = []

from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *

if __name__ == "__main__":

	Load('PandaAnalysisFlat','BTagTreeBuilder')

	builder = root.BTagTreeBuilder()
	builder.SetInput(args.infile)
	builder.SetOutputPath(args.outfile)
	builder.Run()
	builder.Terminate()
