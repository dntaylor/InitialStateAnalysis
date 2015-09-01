InitialStateAnalysis Analyzers
==============================

Running analyzers
-----------------

Analyzers can be run from the [run.py](./scripts/run.py) script. Alternatively one can run over a single analyzer
here.

Creating new analyzers
----------------------

New analyzers should be created from the [AnalyzerTemplate.py](./python/AnalyzerTemplate.py):

```
cp AnalyzerTemplate.py AnalyzerNEW.py
```

From here, you must add the following necessary objects and functions:

* `self.channel`
* `self.final_states`
* `self.initial_states`
* `self.other_states` (optional)
* `self.object_definitions`
* `self.cutflow_labels` (optional)
* `choose_objects(self,rtrow)`
* `choose_alternative_objects(self,rtrow,state)` (optional)
* `preselection(self,rtrow)`
* `selection(self,rtrow)`
*  getIdArgs(self,type)

Look at [AnalyzerWZ.py](./python/AnalyzerWZ.py) for an example.

After adding a new Analyzer, be certain to update [run.py](./scripts/run.py) to make the new analyzer accessible.

Ntuple Content
--------------

The output ntuple content is partially defined by your `initial_states` variable (and optionally `the other_states`
variable). The common ntuple elements are:

`ntuple.root`
* `cutflow` (`TH1`)
* Analysis label (`TTree`, e.g. `WZ`)
  * `select` (`TBranch`)
  * `event` (`TBranch`)
  * `channel` (`Char_t[8]`)
  * `finalstate` (`TBranch`)
  * Final state objects (`TBranch`, e.g. `l1`)
  * Final state object flavors (`Char_t[2]`, e.g. `l1Flv`)
  * ...
  * User defined objects (`TBranch`, e.g. `w`)
  * User defined object flavors (`Char_t[3]`, e.g. `wFlv`)
  * ...

Each of these branches are defined below:

`select`
* `passTight`
* `passLoose`
* `passNM` (NM = numbers 0 to number of objects in final state 1 for each object type)

`event`
* `evt`
* `run`
* `lumi`
* `nvtx`
* `lep_scale`
* `pu_weight`

`finalstate`
* `mass`
* `met`
* `metPhi`
* `jetVetoN` (N = 20,30,40)
* `muonVetoN` (N = 5,10Loose,15)
* `elecVetoN` (N = 10)

Final state object
* `Pt`
* `Eta`
* `Phi`
* `Iso`
* `Chg`

User defined object
* `mass`
* `sT`
* `dPhi`
* `PtN` (N = 1,...,number of objects in composite)
* `EtaN`
* `PhiN`
* `ChgN`
* `met` (If object contains MET)
* `metPhi` (If object contains MET)
