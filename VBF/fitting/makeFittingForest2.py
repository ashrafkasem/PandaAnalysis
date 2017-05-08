#!/usr/bin/env python
from re import sub
from sys import argv,exit
from os import path,getenv
from glob import glob
import argparse
from pprint import pprint

parser = argparse.ArgumentParser(description='make forest')
parser.add_argument('--region',metavar='region',type=str,default=None)
args = parser.parse_args()
out_region = args.region
region = out_region.split('_')[0]
if region=='test':
    is_test = True 
    region = 'signal'
else:
    is_test = False
sname = argv[0]

argv=[]
import PandaAnalysis.Flat.fitting_forest as forest 
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions # kinematics
import PandaAnalysis.Monotop.PandaSelection as sel

basedir = getenv('PANDA_FLATDIR')+'/'
lumi = 35900

def f(x):
    return basedir + x + '.root'

vmap = {}
mc_vmap = {'genBosonPt':'genBosonPt'}
if region in ['signal','test']:
    u,uphi, = ('pfmet','pfmetphi')
elif 'photon' in region:
    u,uphi = ('pfUAmag','pfUAphi')
elif 'single' in region:
    u,uphi = ('pfUWmag','pfUWphi')
elif 'di' in region:
    u,uphi = ('pfUZmag','pfUZphi')
vmap['met'] = u 
vmap['mjj'] = 'jot12Mass'
vmap['deta'] = 'jot12DEta'
vmap['dphi'] = 'jot12DPhi'
weights = {'nominal' : sel.weights[region]%lumi}


factory = forest.RegionFactory(name = region if not(is_test) else 'test',
                               cut = sel.cuts[region],
                               variables = vmap, 
                               mc_variables = mc_vmap, 
                               mc_weights = weights)



if is_test:
    factory.add_process(f('Diboson'),'Diboson')
else:
    if 'signal' in region:
        factory.add_process(f('ZtoNuNu'),'Zvv')
        factory.add_process(f('ZtoNuNu_EWK'),'ewkZvv')
    else:
        if region=='dimuon':
            factory.add_process(f('ZJets'),'Zll',extra_weight)
        else:
            factory.add_process(f('ZJets'),'Zll')
        factory.add_process(f('ZJets_EWK'),'ewkZll')
    factory.add_process(f('WJets'),'Wlv')
    factory.add_process(f('WJets_EWK'),'ewkWlv')
    factory.add_process(f('TTbar'),'ttbar')
    factory.add_process(f('SingleTop'),'ST')
    factory.add_process(f('Diboson'),'Diboson')
    factory.add_process(f('QCD'),'QCD')
    if 'electron' in region:
        factory.add_process(f('SingleElectron'),'Data',is_data=True,extra_cut=sel.triggers['ele'])
    else:
        factory.add_process(f('MET'),'Data',is_data=True,extra_cut=sel.triggers['met'])
    if 'signal' in region:
        factor.add_process(f('ggFHinv_m125'),

factory.run(basedir+'/fitting/fittingForest_%s.root'%out_region)

