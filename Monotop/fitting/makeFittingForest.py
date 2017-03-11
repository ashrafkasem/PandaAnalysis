#!/usr/bin/env python
from re import sub
from sys import argv,exit
from os import path,getenv
from glob import glob
import argparse
parser = argparse.ArgumentParser(description='make forest')
parser.add_argument('--region',metavar='region',type=str,default=None)
region = parser.parse_args().region
if region=='test':
    is_test = True 
    region = 'signal'
else:
    is_test = False

argv=[]
import PandaAnalysis.Flat.fitting_forest as forest 
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions # kinematics
import PandaAnalysis.Monotop.CombinedSelection as sel

basedir = getenv('PANDA_FLATDIR')+'/'
lumi = 35900

def f(x):
    return basedir + x + '.root'

def shift_btags():
    shifted_weights = {}
    if not any([x in region for x in ['signal','top','w']]):
        return shifted_weights 
    for shift in ['BUp','BDown','MUp','MDown']:
        for cent in ['sf_btag','sf_sjbtag']:
            shiftedlabel = '_'
            if 'sj' in cent:
                shiftedlabel += 'sj'
            if 'B' in shift:
                shiftedlabel += 'btag'
            else:
                shiftedlabel += 'mistag'
            if 'Up' in shift:
                shiftedlabel += 'Up'
            else:
                shiftedlabel += 'Down'
            weight = sel.weights[region+'_'+cent+shift]%lumi
            shifted_weights[shiftedlabel] = weight
    return shifted_weights


vmap = {'top_ecf_bdt':'top_ecf_bdt'}
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

weights = {'nominal' : sel.weights[region]%lumi}
weights.update(shift_btags())

factory = forest.RegionFactory(name = region if not(is_test) else 'test',
                               cut = sel.cuts[region],
                               variables = vmap, 
                               mc_variables = mc_vmap, 
                               mc_weights = weights)



if is_test:
    factory.add_process(f('Diboson'),'Diboson')
elif region=='photon':
    factory.add_process(f('GJets'),'Pho')
    factory.add_process(f('SinglePhoton'),'Data',is_data=True)
    factory.add_process(f('SinglePhoton'),'QCD',is_data=True,extra_weights='sf_phoPurity')
else:
    factory.add_process(f('ZtoNuNu'),'Zvv')
    factory.add_process(f('ZJets'),'Zll')
    factory.add_process(f('WJets'),'Wlv')
    factory.add_process(f('TTbar'),'ttbar')
    factory.add_process(f('SingleTop'),'ST')
    factory.add_process(f('Diboson'),'Diboson')
    factory.add_process(f('QCD'),'QCD')
    if 'electron' in region:
        factory.add_process(f('SingleElectron'),'Data',is_data=True)
    else:
        factory.add_process(f('MET'),'Data',is_data=True)
    if region=='signal':
        signal_files = glob(basedir+'/Vector*root')
        for f in signal_files:
            fname = f.split('/')[-1].replace('.root','')
            signame = fname
            replacements = {
                'Vector_MonoTop_NLO_Mphi-':'',
                '_gSM-0p25_gDM-1p0_13TeV-madgraph':'',
                '_Mchi-':'_',
            }
            for k,v in replacements.iteritems():
                signame = signame.replace(k,v)
            factory.add_process(f,signame)


if is_test:
    region = 'test'
factory.run(basedir+'/fitting/fittingForest_%s.root'%region)

