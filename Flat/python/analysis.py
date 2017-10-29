#!/usr/bin/env python

from PandaCore.Tools.Load import Load 
from PandaCore.Tools.Misc import PInfo, PError
import ROOT as root

Load('PandaAnalyzer')

def _analysis(name, verbose, **kwargs):
    if verbose:
        PInfo('PandaAnalysis.Flat.analysis','Summary of analysis %s:'%(name))
    a = root.Analysis(name)
    for k,v in kwargs.iteritems():
        if not hasattr(a, k):
            PError('PandaAnalysis.Flat.analysis','Could not set property %s'%k)
            return None 
        if verbose:
            PInfo('PandaAnalysis.Flat.analysis','    %20s : %s'%(k, 'True' if bool(v) else 'False'))
        setattr(a, k, bool(v))
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
