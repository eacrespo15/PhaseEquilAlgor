# Phase equilibrium & Thermophysical Properties

## About the author:

Emanuel A. Crespo\
PhD in Chemical Engineering: 2017-2021\
University of Aveiro/CICECO Aveiro Institute of Materials

## Introduction:

This repository contains some of the phase equilibrium algorithms (e.g. bubble point calculations, solubility, flash routines, etc) and calculation methods for thermophysical properties (e.g. viscosity calculations) that I have used for some of my PhD works, when using advanced molecular based EoSs, derived from the Statistical Associating Fluid Theory, for the thermodynamic modelling of different systems.

I have used some of these algorithms with a matlab implementation of the soft-SAFT EoS proposed by [Blas and Vega](https://doi.org/10.1080/002689797170707) in 1997. I can not share the code for the EoS itself, as my MATLAB implementation was based on the proprietary FORTRAN code developed by Lourdes F. Vega and co-workers. Nevertheless I can try to help you if you are trying to implement a SAFT-type EoS yourself and want to use some of these codes.

Most of these routines call, at some stage, a routine named `simplex_fug_saft` which was the routine where I calculate the natural logarithm of the fugacity coefficients and the density of my system.

The routines present in this repository are well commented both at the top of the file for input and ouput information, and throughout the code for more detailed information.

## Algorithms in the repository:

### VLE_PTFLASH.m

This file contains a two-phase PT FLASH that can be used with or without acceleration methods either through the Dominant Eigen Value Method (DEVM) or through extrapolation steps. Calls both a Rachford-Rice equations solver (`RR_solver.m`) and a stability analysis for the system (`stability_analysis.m`).

### bubble_point.m

This routine contains an algorithm for bubble point calculations (both bubble point pressures and bubble temperatures) without the need for explicit derivatives of the natural logarithm of the fugacity coefficients.

### MF_FLASH_GUPTA.m

Multiphase flash algorithm based on the work of Gupta et al., Fluid Phase Equilibria, 63 (1991), 65-89. This function calls the `obj_MF_Flash_GUPTA.m` that contains the objective function used by the algorithm.

### MF_FLASH_HEIDEMANN.m

Multiphase flash algorithm with stability analysis based on the works of R. Heidemann from University of Calgary. This function calls `chkphs.m` for stability analysis of the system and `fracts.m` for the material balances of the algorithm.

### MF_FLASH_MICHELSEN.m

Multiphase flash algorithm based on the amazing book of Michelsen & Mollerup - "Thermodynamic models - Fundamentals & Computational Aspects" which I can not recommend enough for those interested in the efficient implementation of thermodynamic models.

Unfortunately, this algorithm still has to be optimized and is still lacking the stability analysis. This routine calls a multiphase Rachford-Rice equations solver in `MPh_Rachford_Rice.m`.

### LLE_Binary and LLE_Ternary.m

Simpler routines for the determination of the Liquid-Liquid Equilibrium (LLE) of binary and ternary systems using the isofugacity criterion.

### SLE_binary.m

Routine for solubility calculations in a binary system. It includes the option to include the effect of heat capacity difference between the solid and liquid phases or the existence of solid-solid transitions. To extend for multicomponent systems in the future.

### SSViscosity.m

This routine implements the Free Volume Theory that I have used coupled to the soft-SAFT EoS for viscosity calculations. As reference, check the work of Llovell et al (J. Phys. Chem. B 2013, 117,8159-8171) [here](https://pubs.acs.org/doi/abs/10.1021/jp401307t). This implementation allows the use of three different mixing rules for the dense term (check the comments throughout the routine for more information).
