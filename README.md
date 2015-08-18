InitialStateAnalysis
======================

The Initial State Analysis (ISA) framework uses ntuples produced with the 
<a href="https://github.com/uwcms/FinalStateAnalysis">FinalStateAnalysis</a> framework.
This framework constructs user defined initial states from final state objects and stores
interesting variables in an output ntuple for further selection.

Installation
------------

ISA requires HiggsCombine tool to produce limits and ROOT.
```
cmsrel CMSSW_7_4_9
cd CMSSW_7_4_9/src
cmsenv
git cms-init
cd recipe
./recipe.sh
cd $CMSSW_BASE/src
scram b -j 16
```

Analyzing data
--------------
The primary analyzer is accessed via the [run.py](Analyzers/scripts/run.py) command. This command has paths to ntuples stored
for convenient access. For example, to run the TT channel of the WZ analysis over all MC samples:

```
# Usage: run.py [analysis] [channel] [period] samples (unix wildcards allowed)
run.py WZ WZ 13 W* T* DY* Z*
```

Jobs can be submitted to the cluster using the `--submit` option:

```
./run.py --submit --jobName=testSubmit Hpp3l Hpp3l 13 D* T* W* Z* 
```

Plotting
--------

Plotting can be accomplished via the [mkplots.py](Plotters/scripts/mkplots.py) command:

```
# Usage: mkplots.py [analysis] [channel] [period] [options]
mkplots.py Hpp3l Hpp3l 13
```

Limits
------

Limits can be run via [mklimits.py](Limits/scripts/mklimits.py). This produces datacards able to be read by the 
`HiggsAnalysis/CombinedLimit` module.

The [mklimits.py](Limits/scripts/mklimits.py) script can produce limits using three different methods: a purely MC driven
method that estimates background from MC samples, a data-driven method with a user defined
sideband and signal region, and a fakerate method (requires the fakerate option on the ntuple production, TODO).

```
# Usage: mklimits [analysis] [region] [period] [options]
mklimits.py Hpp3l Hpp3l 13
```

The datacards can then be processed with the [processdatacards.py](Limits/scripts/processdatacards.py) script:

```
# Usage: processdatacards.py [analysis] [region] [period] [options]
processdatacards.py Hpp3l Hpp3l 13
```

And finally, the limits can be plotted with [plotlimits.py](Plotters/scripts/plotlimits.py):

```
# Usage: plotlimits.py [analysis] [region] [period] [options]
plotlimits.py Hpp3l Hpp3l 13 -bp ee100
```
