# Private signal production

## Setup

First, set up a full installation of `PandaProd`. `PandaUtilities` and `PandaCore` should already be installed. Make sure to call `PandaProd/setuprel.sh` so that everything gets properly checked out. For example:

```bash
cmsrel CMSSW_8_0_20
cd CMSSW_8_0_20/src
cmsenv
git cms-init
git clone https://github.com/PandaPhysics/PandaProd
git clone https://github.com/PandaPhysics/PandaUtilities
git clone https://github.com/PandaPhysics/PandaCore
sh PandaProd/setuprel.sh
git clone https://github.com/sidnarayanan/PandaAnalysis
scram b -j16
```

## Production

The configuration used is `PandaProd/Ntupler/test/runNtuplerCondor.py`. This is identical to `runNtupler.py`, but picks up input/output files differently.

The main submission code is `PandaAnalysis/Monotop/prod`. To run, you will need a catalog of files on the T2. You can find them here [1]. Copy `cfg/` into `PandaAnalysis/Monotop/prod`. Then, to submit, simply call `python submit.py <model name>`. For example, you can do:

```bash
ls cfg | sed "s?.txt??g" | xargs -n 1 python submit.py
```
 
The output files will land in `~/cms/hist/monotop_private_panda`. Note, there is not really a good way of resubmitting failed jobs right now. 

[1] `/home/snarayan/cms/cmssw/analysis/CMSSW_8_0_20/src/PandaAnalysis/Monotop/prod/cfg` 

## Catalog

Like with regular MC production, you need to produce a catalog of the output. There is an analogous script in `config/catalogProd.py`. You can run it:

```bash
./catalogProd.py --outfile ~/public_html/histcatalog/$(date +%Y%m%d).cfg
```

## Submit analyzer

Running the analyzer is the same as before, just make sure to update `T3/setup.sh` so that the catalog points to the private production catalog. The skimming configuration script can automatically differentiate between remote and local data, so nothing else needs to change in the skimming and merging steps.

