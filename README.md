![alt text](common/llogo.png)

_A tool to calculate thermodynamic contributions from ORCA with full control_

***
## Introduction

**otherm** provides the functionality to calculate entropic (_S_), enthalpic (_H_) and Gibbs free (_G_) energies from
ORCA output files with control beyond that provided in ORCA. Different standard states (1 atm / 1 M) can be specified,
symmetry numbers calculated, free energy calculations at different temperatures from the same output performed, and 
alternate methods used for calculating the entropic contribution used.


## Output

```bash
\----------------------------------------------------------------------------------
Filename                                                                   test.out
Temperature (K)                                                               298.0
Standard state is                                                               1 M
Calculating using the method of                                              Grimme
                                                       Chem. Eur. J. 2012, 18, 9955

Symmetry number (σ)                                                               1
Molecular weight (amu)                                                       136.15

Total entropy (J K-1 mol-1)                                                  356.90
Total enthalpy (J mol-1)                                             -1205452153.88
Total free energy (J mol-1)                                          -1205558509.65

For convenience E, H, G in Hartrees                                                
-459.286861041774,-459.13247760703814,-459.1729863786075
----------------------------------------------------------------------------------
```


## Usage

See _examples_/ for further example usage and

```bash
python otherm.py --help
```

for help.

## Dependencies
* Python 3
* numpy