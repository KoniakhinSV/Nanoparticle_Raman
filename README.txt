Raman Spectra Analysis and Simulation Tool for Nanodiamonds
============================================================

This software provides Python-based tools for processing experimental Raman spectra of diamond nanoparticles and calculating theoretical spectra from user-defined nanodiamond powder parameters.

It combines physical modeling approaches with practical fitting spectra with reconstruction of nanodiamond size histogram
   - raman_fit_NN.py
   - raman_fit_Metropolis.py
and simulation capabilities
   - raman_calculate.py

------------------------------------------------------------------
Physical Models Behind the Algorithms
------------------------------------------------------------------

The tool implements recent theoretical developments in Raman spectroscopy of nanodiamonds, particularly:

1. Discrete size-quantized phonon modes in finite-size crystallites
   Due to finite crystallite size, standing wave like phonons (non-zero effetive wavevectors) contribute to Raman scattering.
   This results in a downshift and asymmetric broadening of the diamond peak compared to bulk diamond due to negative mass dispersion.

2. Continuous approach: "elasticity theory" for optical phonons
   Discrete (atomistic) dynamical matrix method can be mapped to continuous Klein-Fock-Gordon equation with Euclidean metric (scalar "elasticity theory" for optical phonons). Raman cross section within bond polarization model can be mapped to the square of the phonon wave-function integral over nanoparticle volume (form-factor).

3. Lattice Disorder Effects
   Point defects and impurities, such as vacancies or substitutional atoms, cause additional shifts
   and broadening of the Raman line. For nanodiamonds, these contributions are of the same order
   of magnitude as size-confinement effects.

Accurate spectral analysis requires considering *both* finite-size confinement and impurity-related contributions.

------------------------------------------------------------------
Features
------------------------------------------------------------------

1. Experimental Raman Spectrum Fitting

   - Neural Network Method (raman_fit_NN.py)
     Uses a pre-trained neural network to predict nanoparticle size distribution and broadening parameters. Works faster.

   - Metropolis Optimization (raman_fit_Metropolis.py)
     Iteratively refines the size distribution to minimize mean squared deviation between experimental and simulated spectra. Works better for multi-peak size ditributions.

   Both methods:
     - Validate input spectra
     - Apply the combined physical model (confinement + impurity effects)
     - Output a ZIP archive containing reconstructed spectra, histograms, and estimated parameters (C_imp, Gamma0)

2. Spectrum Calculation from diamond nanopowder parameters*

   - Forward spectrum modeling (raman_calculate.py)
     Generates Raman spectra from:
       D      : Mean particle diameter (nm)
       dD     : Size dispersion (log-normal distribution assumed)
       Cimp   : Impurity concentration
       Gamma0 : Intrinsic linewidth
---------
  * At their own risk, the User can change the dispersion paramters A,B in the "raman_spectrum_single_size" function from raman_routines module.

------------------------------------------------------------------
Input Data Formats
------------------------------------------------------------------

### 1. Neural Network and 2. Metropolis spectrum fitting:
  input.csv with two columns (no header):
    Wavenumber (cm^-1), Intensity

  Requirements:
    - Wavenumber range covers at least 1300–1340 cm^-1 with >= 40 points in this interval.
    - Wavenumber values should be evenly destributed (no too small and too large intervals).
    - Wavenumber values strictly increasing.
    - Intensities are normalized automatically.

    Example:
1299.6,0.111
1300.0,0.123
1300.4,0.134
...
1340.0,0.115
1340.4,0.0

### 3. Forward spectrum calculation:
  parameters.csv (header line and parameter line):
    D,dD,Cimp,Gamma0
    4.0,0.3,1.2,2.5

  Constraints:
    1.5 <= D <= 9.0
    0.05 <= dD <= 1.0
    0.0 <= Cimp <= 5.0
    0.5 <= Gamma0 <= 15.0


Examples of input.csv and parameter.csv files are presented in the Repository.

------------------------------------------------------------------
Usage
------------------------------------------------------------------

### 1. Neural Network Fitting
1. Place your experimental spectrum in the same directory as raman_fit_NN.py and name it input.csv.
2. Run:
   python raman_fit_NN.py
   (Run.sh file gives an example of launch in Linux with conda environment)

3. The script will produce:
     - result_Raman_Fit_NN.zip containing:
     - NN_spectra.csv – Experimental vs fitted spectrum
     - NN_histogram.csv – Reconstructed particle size distribution
     - NN_spectra.png – Plot of spectrum match
     - NN_histogram.png – Plot of size distribution

### 2. Metropolis fitting:
   Works the same as Neural Network-based approach. Produces result_Raman_Fit_Metropolis.zip

### 3. Forward calculation:
1. Tune the diamond nanopowder paratmertes in the parameters.csv file of the following format:
  D,dD,Cimp,Gamma0
  4.0,0.3,1.2,2.5
2. Run:
   python raman_calculate.py
3. The script will produce:
     - result_Raman_Calculated.zip containing
     - NN_spectra.csv – Experimental vs fitted spectrum
     - NN_histogram.csv – Reconstructed particle size distribution
     - NN_spectra.png – Plot of spectrum match
     - NN_histogram.png – Plot of size distribution

------------------------------------------------------------------
Output ZIP Contents
------------------------------------------------------------------

Each ZIP archive includes:
  *_spectra.csv   - Experimental vs. reconstructed/calculated spectra
  *_spectra.png   - Plot of spectra comparison
  *_hist.csv      - Predicted or input size distribution
  *_hist.png      - Histogram of particle sizes

------------------------------------------------------------------
Dependencies
------------------------------------------------------------------
- Python >= 3.8
- numpy
- matplotlib

Install with:
  pip install numpy matplotlib

------------------------------------------------------------------
Repository Structure
------------------------------------------------------------------
raman_routines.py         - Core spectrum modeling, fitting, and validation functions
raman_fit_NN.py           - Neural Network-based fitting
raman_fit_Metropolis.py   - Metropolis-based fitting
raman_calculate.py        - Forward spectrum calculation
model_weights_0.npz       - Pre-trained NN weights (required for NN fitting)

------------------------------------------------------------------
License
------------------------------------------------------------------
Released under the MIT License.
