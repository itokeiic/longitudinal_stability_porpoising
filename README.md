# longitudinal_stability_porpoising
This repository contains codes used in the journal article titled ["Longitudinal Stability Augmentation of Seaplanes in Planing"](https://www.researchgate.net/publication/303504143_Longitudinal_Stability_Augmentation_of_Seaplanes_in_Planing)

The main program is suspended_planing.py.  It will output plots of maximum real part of eigen values of plant matrices of rigidly-supported and flexibly-supported planing crafts (seaplanes) vs. Froude numbers (non-dimensional planing speed). It essentially tells you at what speed the crafts become unstable (i.e. the real eigen values become positive).  Different figures corresponds to different spring stiffness (k_f for front support, k_b for back support) and damping coefficients (c_f for front support, c_b for back support).  It will also produce a figure of plots showing trim angles vs. Froude number with respect to two center of gravity location in rigidly-supported craft. Please refer the paper above for details.

It is written in Python 2.7 and can be run at the command prompt as

```
python suspended_planing.py
```

The program requires numpy, scipy, matplotlib, and LaTeX installation.