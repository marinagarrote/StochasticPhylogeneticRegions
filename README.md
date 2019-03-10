# StochasticPhylogeneticVarieties

This repository contains the code to compute the distance to the stochastic phylogenetic varieties, for the case presented in the paper *Computing the distance to the positive part of phylogenetic varieties* by M. Casanellas, J. Fernández-Sánchez and M. Garrote López submitted as a Computation Presentation at the **MEGA 2019 conference**.

## Prerequisites

[**Macaulay2**](http://www.math.uiuc.edu/Macaulay2/) and [**PhCpack**](http://homepages.math.uic.edu/~jan/phcpy_doc_html/welcome.html) should be previously installed.

## Running the tests

File *AlgorithmStochasticPhylogeneticVarieties.m2* contains the main code.

### Saturation ideal

The following lines are comented in the code since the execution of them could take several days. Uncomment them in case you want to compute the degree for different parameters.

```
J = saturate(ideal DFu, ideal(x_0*x_1*x_2*x_3*x_4));
degJ = degree J;
```

### Test the algorithm for parameters k0, m0

Substitute parameters k0 and m0 for desired values

```
k0 = 2/10;
m0 = 14/10;
```


