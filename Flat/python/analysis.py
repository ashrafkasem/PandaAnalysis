#!/usr/bin/env python

from PandaCore.Tools.Load import Load 
from PandaCore.Tools.Misc import PInfo, PError
import ROOT as root

Load('PandaAnalyzer')

def _dump(a):
    PInfo('PandaAnalysis.Flat.analysis','Summary of analysis %s:'%(a.name))
    for k in dir(a):
        if k[0] == '_':
            continue
        PInfo('PandaAnalysis.Flat.analysis','    %20s = %s'%(k, 'True' if bool(getattr(a, k)) else 'False'))



def _analysis(name, verbose, **kwargs):
    a = root.Analysis(name)
    for k,v in kwargs.iteritems():
        if not hasattr(a, k):
            PError('PandaAnalysis.Flat.analysis','Could not set property %s'%k)
            return None 
        setattr(a, k, bool(v))
    setattr(a, 'dump', lambda : _dump(a))
    if verbose:
        a.dump()
    return a

def analysis(name, **kwargs):
    return _analysis(name, verbose=True, **kwargs)


# predefined!
monotop = lambda v=False : _analysis(
        name = 'monotop',
        verbose = v,
    )

vbf = lambda v=False : _analysis(
        name = 'vbf',
        verbose = v,
        vbf = True,
        fatjet = False,
        btagSFs = False,
        puppi_jets = False
    )

monoh = lambda v=False : _analysis(
        name = 'monoh',
        verbose = v,
        monoh = True,
    )

gghbb = lambda v=False : _analysis(
        name = 'gghbb',
        verbose = v,
        monoh = True,
        recoil = False,
        ak8 = True,
    )
wlnhbb = lambda v=False : _analysis(
        name = 'wlnhbb',
        verbose = v,
        monoh = True,
        recoil = False,
        ak8 = False,
        fatjet = False,
        btagSFs = False,
        btagWeights = True,
        useCMVA = True,
        complicatedLeptons = True,
        genPartonStudy = True
    )
