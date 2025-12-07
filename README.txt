============================================================
Raman Spectra Analysis and Simulation Tool for Nanodiamonds
============================================================



==================================================================
QUICK START
==================================================================

Available both as a script repository and an online tool:
1. https://github.com/KoniakhinSV/Nanoparticle_Raman
2. https://nanoraman.pythonanywhere.com

This toolkit provides Python-based methods for:

- Fitting experimental Raman spectra of nanodiamonds
- Simulating Raman spectra from user-defined nanodiamond parameters

Two fitting approaches are available:
   • raman_fit_NN.py — neural-network-based
   • raman_fit_Metropolis.py — Metropolis optimization

Simulation is performed by:
   • raman_calculate.py

Theory: https://arxiv.org/abs/2511.23053


------------------------------------------------------------------
Usage
------------------------------------------------------------------

### 1. Neural Network Fitting
1. Place input spectrum as **input.csv** next to raman_fit_NN.py
2. Adjust **parameters.csv** if needed
3. Run:
   `python raman_fit_NN.py`
4. Output: **result_Raman_Fit_NN.zip**

### 2. Metropolis Fitting
1. Place input spectrum as **input.csv** next to raman_fit_Metropolis.py
2. Edit **parameters.csv** for initial values and toggles
3. Run:
   `python raman_fit_Metropolis.py`
4. Output: **result_Raman_Fit_Metropolis.zip**

### 3. Forward Spectrum Calculation
1. Set D, dD, Cimp, Gamma0 in **parameters.csv**
2. Run:
   `python raman_calculate.py`
3. Output: **result_Raman_Calculated.zip**



==================================================================
DESCRIPTION
==================================================================

------------------------------------------------------------------
Physical Models Behind the Algorithms
------------------------------------------------------------------

The implementation follows modern theories of Raman scattering in nanodiamonds. Accurate reconstruction requires both confinement and impurity effects:

1. **Finite-size phonon confinement**
   Discrete quantized phonon modes cause peak downshift and asymmetric broadening
   (arXiv:1803.01653v2; https://doi.org/10.1021/acs.jpcc.8b05415)

2. **Continuum “elasticity-like” optical phonon model**
   Atomistic dynamics mapped to a Klein–Fock–Gordon–type equation
   (arXiv:1806.08100v1; https://doi.org/10.1021/acs.jpcc.8b07061)

3. **Lattice disorder and impurities**
   Vacancies/substitutional atoms introduce size-independent shifts and broadening
   (arXiv:2403.17310; https://doi.org/10.1016/j.diamond.2024.111182)

4. **Raman Spectra Fitting procedure**
Neural Network and Metropolis implementation
(https://arxiv.org/abs/2511.23053)


------------------------------------------------------------------
Features
------------------------------------------------------------------

**1. Fitting Experimental Raman Spectra**

- **Neural Network (raman_fit_NN.py)**
  Fast reconstruction of size distribution and impurity parameters.

- **Metropolis Optimization (raman_fit_Metropolis.py)**
  Iterative fitting; better suited for multimodal size distributions.
  Can initialize from parameters.csv.

Both methods:
- Validate input spectra
- Apply the full physical Raman model
- Output a ZIP archive with fitted spectrum, size histogram, and parameters


**2. Forward Spectrum Calculation (raman_calculate.py)**
Generates Raman spectra from nanopowder parameters:
- D (mean size), dD (dispersion), Cimp, Gamma0
(Log-normal distribution assumed; dispersion constants A,B can be edited manually.)


------------------------------------------------------------------
Input Data Formats
------------------------------------------------------------------

------------------------------
parameters.csv
------------------------------
Used by **all tools**.

Fields:
  background    : 1/0 — apply/pass NN-based background subtraction
  D, dD         : Mean particle diameter, size dispersion
  Cimp          : Impurity concentration
  Gamma0        : Intrinsic linewidth
  adjust_Cimp   : 1/0 — allow/prevent Metropolis to vary Cimp
  adjust_Gamma0 : 1/0 — allow/prevent Metropolis to vary Gamma0

Constraints:
  1.5 ≤ D ≤ 9.0
  0.05 ≤ dD ≤ 1.0
  0 ≤ Cimp ≤ 5
  0.5 ≤ Gamma0 ≤ 35

Example:
background,D,dD,Cimp,Gamma0,adjust_Cimp,adjust_Gamma0
1.0,4.0,0.1,0.0,14.5,1.0,1.0


------------------------------
input.csv
------------------------------
Required for **NN** and **Metropolis** fitting.

Two-column file *(no header)*:
  Wavenumber (cm⁻¹), Intensity

Requirements:
- Must fully include 1300–1340 cm⁻¹ with ≥40 points
- Strictly increasing, nearly uniform grid
- Intensities normalized internally

Example:
1299.5,0.111
1300.0,0.123
…
1340.5,0.0


------------------------------------------------------------------
Output ZIP Contents
------------------------------------------------------------------

Each ZIP contains:
  *_spectra.csv   — Experimental vs reconstructed/calculated spectrum
  *_spectra.png   — Spectrum comparison
  *_hist.csv      — Size distribution (fitted or user-defined)
  *_hist.png      — Histogram plot


------------------------------------------------------------------
Dependencies
------------------------------------------------------------------
- Python ≥ 3.9
- numpy
- matplotlib

Install:
`pip install numpy matplotlib`


------------------------------------------------------------------
Repository Structure
------------------------------------------------------------------

EXECUTABLES:
  raman_fit_NN.py
  raman_fit_Metropolis.py
  raman_calculate.py

INPUT:
  parameters.csv
  input.csv

SUPPORT FILES:
  raman_routines.py
  NN_load.npz
  NN_loadB.npz

README.txt


------------------------------------------------------------------
License
------------------------------------------------------------------
MIT License
