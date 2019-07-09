# otherm examples

Provided in this directory are three ORCA output files in which frequency calculation have been performed. The first 
`test.out` is a true minimum energy structure with no imaginary frequencies and the second `test_imag.out` with one 
spurious imaginary frequency and `methane.out` which provides a test for the symmetry evaluation.
 
 
## 1. ORCA vs otherm
To calculate the free energy contribution
using the same method as in ORCA

```bash
python ../otherm.py test.out -ss 1atm -sn 3
```

which will print a total free energy (_G_) of -459.1749664172693 Hartrees, which differs from that calculated in ORCA
(under _Final Gibbs free enthalpy_) by < 0.05 kcal mol-1. 

## 2. Standard states
By default most electronic structure codes by default utilise a 1 atm 'standard state' for the calculation of the
translational partition function. However, for a solution phase reaction the effective volume is dramatically 
overestimated. To partially alleviate this deficiency a 1 M 'standard state' was proposed in with a scaled down
translational component (e.g. see [here](https://doi.org/10.1021/jp205508z)). To calculate the free energy at a 1M
standard state

```bash
python ../otherm.py test.out -ss 1M
```
or 
```bash
python ../otherm.py test.out
```
as by default **otherm** will calculate the free energy in a 1M standard state


## 3. Symmetry numbers
In ORCA 4.0 a symmetry number (σ) of 3 was assumed for all molecules, which is obviously not correct for a molecule with 
_C1_ symmetry. Instead *otherm* will calculate the symmetry number based on the number of unique rotation axes e.g.
for methane we have σ=12 (see [here](https://doi.org/10.1007/s00214-007-0328-0))

```bash
python ../otherm.py methane.out
```

will print find σ=12.

Due to the number of rotation axes that are searched in the algorithm the symmetry number calculation can become rather
slow for molecules with more than 50 atoms. Therefore, by default σ=1 for molecules with >50 atoms unless the 
`--calc_sym` flag is given. Alternatively, the symmetry number can be specified manually with `-sn` e.g.

```bash
python ../otherm.py methane.out -sn 12
```

will produce an identical output to the above.

## 4. Temperature
By default the free energy contribution in ORCA is calculated at 298 K, requiring a second frequency calculation if
a new free energy at a different temperature is required. In *otherm* a new temperature is invoked with `-t` or
`--temp` and an temperature in K i.e.

```bash
python ../otherm.py test.out -t 373
```

## 5. Different Methods
In ORCA the default method for the entropy calculation is that of [Grimme](https://doi.org/10.1002/chem.201200497) the 
quasi rigid rotor harmonic oscillator (qRRHO) method, where the low frequency normal modes are treated as rotors as is
the default in *otherm*. To use the ideal gas model (IGM) as in Gaussian

```bash
python ../otherm.py test.out -m igm
```

or the method of [Truhlar](https://doi.org/10.1021/jp205508z) where the low frequency modes are scaled to a constant
value defined by `-s` (default to 100 cm-1)


```bash
python ../otherm.py test.out -m truhlar
```

To use Grimme's method with non-default parameters (alpha, w_0 in the [paper](https://doi.org/10.1002/chem.201200497))

```bash
python ../otherm.py test.out -m grimme --alpha 4.5 --w0 50
```
where the defaults are : alpha = 4.0 and w0 = 100 cm-1.

