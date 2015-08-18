Limits
======

Producing datacards
-------------------

Datacards are produced via the [mklimits.py](./scripts/mklimits.py) script. The produced
datacards are then able to be processed via the Higgs Combine tool. The structure of
the datacards is defined below. The [datacard.py](./python/datacard.py) class stores the information
about the signal and backgrounds and the [Limits.py](./python/Limits.py) class creates the datacard
from information (such as signal branching fractions, full selections, background estimation
method) from the [mklimits.py](./scripts/mklimits.py) script.

Datacard structure
------------------

A sample datacard is copied below:

```
#Hpp3l
imax  1 number of channels
jmax  6 number of backgrounds
kmax  7 number of nuisance parameters
------------
bin 1
observation 0
------------
bin                   1          1          1          1          1          1          1
process            hpp500       ttz        top        ttw        wz         zz         tt
process               0          1          2          3          4          5          6
rate              8.005e+00  1.396e-02  1.000e-10  1.000e-10  1.000e-10  1.000e-10  1.000e-10
------------
lumi        lnN     1.026      1.026      1.026      1.026      1.026      1.026      1.026
eid         lnN     1.01       1.01       1.01       1.01       1.01       1.01       1.01
muid        lnN     1.005      1.005      1.005      1.005      1.005      1.005      1.005
muiso       lnN     1.002      1.002      1.002      1.002      1.002      1.002      1.002
sigmc       lnN     1.15         -          -          -          -          -          -
mcdata      lnN       -          -          -          -        1.056      1.105      1.024
pdf         lnN       -        1.105        -        1.289        -          -          -
```
