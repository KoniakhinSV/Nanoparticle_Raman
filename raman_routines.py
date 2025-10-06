import numpy as np
import matplotlib
matplotlib.use("Agg") # uncomment for standalone
import matplotlib.pyplot as plt

import os
import zipfile
import copy

#===================================================================

def find_load_npz_file_alt():

    #cwd = os.getcwd()
    cwd = os.path.dirname(os.path.abspath(__file__))
    
    for filename in os.listdir(cwd):
        if filename.endswith('_load.npz'):
            return filename


#===================================================================

#===================================================================
def define_constans():
# Constants
    OMEGA_MIN, OMEGA_MAX, OMEGA_POINTS = 1300, 1340, 100
    OMEGA_AXIS = np.linspace(OMEGA_MIN, OMEGA_MAX, OMEGA_POINTS)

    # Fixed parameters
    C_IMP = 1.2*2
    GAMMA0 = 2.5*2

    BINS = np.arange(1.5, 15.01, 0.25)
    #BINS = np.arange(1.5, 14.01, 0.5)

    N_BINS = len(BINS)
    CORRECTION_COEFF = 0.05

    my_dict = {"NORMALIZE": BINS**3.0,
                "OMEGA_AXIS": OMEGA_AXIS, "BINS": BINS, "N_BINS": N_BINS,
                "CORRECTION_COEFF": CORRECTION_COEFF,"N_Cimp": 16, "N_Gamma0":16}
    return my_dict

CONST = define_constans()

#===================================================================

#===========================================================================

def c_imp_vs_index(ind):
    return 0.2*ind

def Gamma0_vs_index(ind):
    return 1.0 *(1.28**ind)


#===========================================================================

# --- Helper functions from earlier ---

def disorder_effects(c_imp, L_nm):
    delta_omega = -2.6 * c_imp
    Gamma_disorder = 1.5 * c_imp * (3.0 / L_nm)
    return delta_omega, Gamma_disorder

def raman_spectrum_single_size(R_cm, c_imp, Gamma0, n_max=5, CONST = CONST):
    omega_axis = CONST["OMEGA_AXIS"]
    a0 = 0.357e-7
    A, B = 1332.5, 85.0
    #A, B = 1333.0, 66.0
    L_nm = 2 * R_cm * 1e7
    d_omega, Gamma_dis = disorder_effects(c_imp, L_nm)
    total_shift = d_omega
    total_Gamma = Gamma0 + Gamma_dis

    I = np.zeros_like(omega_axis)
    for n in range(1, n_max + 1):
        q = np.pi * n / R_cm
        omega_n = (A + 0) - B * (q ** 2 * a0 ** 2 / 8) + total_shift
        In = 6 * (4 / 3) * np.pi * (R_cm ** 3) / ((np.pi * n) ** 2)
        I += In * (n * total_Gamma / 2) ** 2 / ((omega_axis - omega_n) ** 2 + ( n * total_Gamma / 2) ** 2)
    return I

def generate_lognormal_counts(mu_nm=4.0, sigma_nm=0.3, N_total=1000, CONST = CONST):
    sizes_nm = CONST["BINS"]
    log_sizes = np.log(CONST["BINS"])
    mu = np.log(mu_nm)
    sigma = sigma_nm
    pdf = (1 / (sizes_nm * sigma * np.sqrt(2 * np.pi))) * \
          np.exp(- (log_sizes - mu) ** 2 / (2 * sigma ** 2))
    counts = pdf / pdf.sum() * N_total
    return sizes_nm, counts

def generate_normal_counts(mu_nm=4.0, sigma_nm=0.3, N_total=1000, CONST=None):
    sizes_nm = CONST["BINS"]
    mu = mu_nm
    sigma = sigma_nm
    
    # Normal distribution PDF
    pdf = (1 / (sigma * np.sqrt(2 * np.pi))) * \
          np.exp(- (sizes_nm - mu) ** 2 / (2 * sigma ** 2))
    
    counts = pdf / pdf.sum() * N_total
    return sizes_nm, counts

def total_raman_spectrum(distribution_counts, c_imp, Gamma0, n_max=5, CONST = CONST):
    sizes_nm = CONST["BINS"]
    total_I = np.zeros_like(CONST["OMEGA_AXIS"])

    for size_nm, count in np.column_stack((CONST["BINS"],distribution_counts)):
        R_cm = (size_nm / 2) * 1e-7
        spectrum = raman_spectrum_single_size(R_cm, c_imp, Gamma0, n_max, CONST = CONST)
        total_I += count * spectrum

    return total_I / np.max(total_I)  # normalize

def mean_squared_deviation(spectrum1, spectrum2):
    return np.mean((spectrum1 - spectrum2) ** 2)









#====================================================================================================
#      M E T R O P O L I S


#			M     M
#			MM   MM
#			M M M M
#			M  M  M
#			M     M
#			M     M
#			M     M


#====================================================================================================






def metropolis(known_spectrum, CONST = CONST, PARAMS = False):


    # --- Step 2: Initialize current distribution randomly / take initial params ---
    if PARAMS and isinstance(PARAMS, (list, tuple)) and len(PARAMS) == 4:
        mu0, sigma0, C_IMP, GAMMA0 = PARAMS
    else:
        C_IMP = 0.7
        GAMMA0 = 15.0
        #mu0 = np.random.uniform(2.0, 7.0)
        #sigma0 = np.random.uniform(0.2, 0.6)
        mu0 = 5
        sigma0 = 0.3


    _, current_counts = generate_lognormal_counts(mu_nm=mu0, sigma_nm=sigma0, N_total=1000, CONST = CONST)
    current_counts /= np.max(current_counts)
    current_spectrum = total_raman_spectrum(current_counts, C_IMP, GAMMA0, CONST=CONST)
    current_msd = mean_squared_deviation(current_spectrum, known_spectrum[:,1])

    # --- Step 3: Iterative Metropolis-like refinement ---

    history_msd = [current_msd]
    history_counts = [current_counts.copy()]
    MAX_ITER = 1000 * 3

    for step in range(MAX_ITER):
        # Generate small correction
        mu_c = np.random.uniform(2.0, 7.0)
        #sigma_c = np.random.uniform(0.05 , 0.6)
        #sigma_c = 0.08 * (2.0**np.random.uniform(0.0 , 3.6))
        sigma_c = np.random.uniform(0.05, 0.3*6/mu_c)
        _, correction = generate_lognormal_counts(mu_c, sigma_c, N_total=1000, CONST=CONST)
        correction /= np.max(correction)

        # Add correction
        trial_counts = current_counts + 0.05 * correction * np.random.rand()
        #trial_counts = current_counts + 0.05 * correction * (2*np.random.rand() - 1)
        trial_counts = np.clip(trial_counts,
                                a_min = 0.0, a_max = 40.0);
        trial_counts /= np.max(trial_counts)  # keep normalized
        
        C_IMP_trial = np.clip(C_IMP + np.random.uniform(-0.01,0.01),
                                a_min = 0.0, a_max = 3.0);
        GAMMA0_trial = np.clip(GAMMA0 + np.random.uniform(-0.05,0.05),
                                a_min = 0.5, a_max = 40.0);

        trial_spectrum = total_raman_spectrum(trial_counts, C_IMP_trial, GAMMA0_trial,
        CONST=CONST)
        trial_msd = mean_squared_deviation(trial_spectrum, known_spectrum[:,1])

        delta_msd = trial_msd - current_msd
        rr = np.random.uniform();  #or np.exp(delta_msd/0.01) > rr
        #if np.exp(-delta_msd/0.00001) > rr: uncomment for full Metropolis
        if delta_msd < 0:   # comment for full Metropolis
            current_counts = trial_counts
            current_spectrum = trial_spectrum
            C_IMP = C_IMP_trial
            GAMMA0 = GAMMA0_trial
            current_msd = trial_msd
            history_msd.append(current_msd)
            history_counts.append(current_counts.copy())


    return known_spectrum, np.column_stack((CONST["BINS"],current_counts)), np.column_stack((CONST["OMEGA_AXIS"],current_spectrum)), C_IMP, GAMMA0, history_msd











#=====================================================================

#   F I G U R E S


#   for NN and Metropolis


#			FFFFFFF
#			F
#			F
#			FFFF
#			F
#			F
#			F




#=====================================================================








def create_output_zip(original_csv_path, zip_path, result_main, mode = "NN", CONST = CONST, Cleanup = True):
    # Prepare paths
    known_spectrum, reconstructed_distribution, reconstructed_spectrum, C_IMP, GAMMA0, _ = result_main
    #reconstructed_distribution[:,1] = reconstructed_distribution[:,1]/np.max(reconstructed_distribution[:,1])
    volume_distribution = reconstructed_distribution[:,1]*CONST["NORMALIZE"]
    volume_distribution = volume_distribution/np.max(volume_distribution)

    mean_D = np.sum(reconstructed_distribution[:,0] * reconstructed_distribution[:,1])/np.sum(reconstructed_distribution[:,1])
    mean_SCA = np.sum(reconstructed_distribution[:,0]**4 * reconstructed_distribution[:,1])/np.sum(reconstructed_distribution[:,0]**3 * reconstructed_distribution[:,1])

    base_dir = os.path.dirname(original_csv_path)
    original_csv_name = os.path.basename(original_csv_path)

    spectra_csv_path = os.path.join(base_dir, mode+"_spectra.csv")
    hist_csv_path = os.path.join(base_dir, mode+"_histogram.csv")

    spectra_plot_path = os.path.join(base_dir, mode+"_spectra.png")
    hist_plot_path = os.path.join(base_dir, mode+"_histogram.png")

    # Save processed CSV
    exp_ram = np.column_stack((known_spectrum[:,0], known_spectrum[:,1],reconstructed_spectrum[:,1]))
    np.savetxt(spectra_csv_path, exp_ram, delimiter=",", fmt="%.6f", header='Energy, original, reconstructed')
    
    #exp_hist = np.column_stack((reconstructed_distribution))
    exp_hist = np.column_stack((reconstructed_distribution[:,0],reconstructed_distribution[:,1],volume_distribution))
    np.savetxt(hist_csv_path, exp_hist, delimiter=",", fmt="%.6f", header='Size_nm, By_N, By_Vol')

    #==========================================================================


    # Plot spectrum match
    plt.figure(figsize=(10, 6))
    plt.plot(known_spectrum[:,0], known_spectrum[:,1], label='Known spectrum', linewidth=2, color = 'black')
    plt.plot(known_spectrum[:,0], reconstructed_spectrum[:,1], label='Reconstructed spectrum', linestyle='--',color = 'blue')
    plt.xlabel('ω (cm⁻¹)')
    plt.ylabel('I(ω)')
    plt.title(mode + ': Spectrum reconstruction')
    plt.legend()
    plt.savefig(spectra_plot_path)
    plt.close()

    # Plot distribution match
    plt.figure(figsize=(10, 6))
    #plt.bar(sizes_nm - 0.1, true_counts, width=0.2, label='True distribution')
    #plt.bar(sizes_nm + 0.1, reconstructed_distribution, width=0.2, label='Reconstructed distribution')
    plt.bar(reconstructed_distribution[:,0]-0.05, reconstructed_distribution[:,1], width=0.1, label='By number', color = 'blue')
    plt.bar(reconstructed_distribution[:,0]+0.05, volume_distribution, width=0.1, label='By volume',color='lightgray', edgecolor='blue', linestyle='dotted')
    plt.xlabel('Particle size (nm)')
    plt.ylabel('Normalized counts')
    plt.title(f"Size distribution: $C_{{\\mathrm{{imp}}}} = {C_IMP:.2f}$%  $\Gamma_0 = {GAMMA0:.2f}~\\mathrm{{cm}}^{{-1}}$")
    plt.text(0.0,-0.13,f"$\\langle D \\rangle_{{N}} = ${mean_D:.2f}nm  $\\langle D \\rangle_{{V}} = ${mean_SCA:.2f}nm")

    plt.legend()
    plt.savefig(hist_plot_path)
    plt.close()

    #==========================================================================
    '''# Save plot
    plt.figure()
    plt.plot(processed_array[:, 0], processed_array[:, 1])
    plt.xlabel("X")
    plt.ylabel("Processed Y")
    plt.title("Processed Data Plot")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()'''

    # Create ZIP archive
    with zipfile.ZipFile(zip_path, 'w') as zipf:
        zipf.write(spectra_csv_path, arcname=mode+"_spectra.csv")
        zipf.write(hist_csv_path, arcname=mode+"_hist.csv")

        zipf.write(spectra_plot_path, arcname=mode+"_spectra.png")
        zipf.write(hist_plot_path, arcname=mode+"_hist.png")

    # Cleanup
    if Cleanup:
        os.remove(spectra_csv_path)
        os.remove(hist_csv_path)
        os.remove(spectra_plot_path)
        os.remove(hist_plot_path)






#=====================================================================

#   S A M P L E



#			 SSSSS
#			S     S
#			S
#			 SSSS
#			     S
#			S     S
#			 SSSS




#=====================================================================


def generate_sample(npts = 140, CONST = CONST, n_distributions=3, noiselevel = 0.1, PARAMS = False):
    # --- Step 1: Generate true (unknown) distribution and its spectrum ---

    true_counts = np.zeros(CONST["N_BINS"])

    # Combine several random log-normal distributions

    if PARAMS == False:
        for _ in range(n_distributions):
            mu = np.random.uniform(2.0, 6.0)
            #sigma = np.random.uniform(0.2 - 0.12, 0.6 - 0.3)
            sigma = np.random.uniform(0.05, 0.3*6/mu)
            _, counts = generate_lognormal_counts(mu_nm=mu,
            sigma_nm=sigma, N_total=np.random.uniform(300, 800), CONST = CONST)
            true_counts += counts
        C_IMP = np.random.uniform(0.0,2.7)
        GAMMA0 = 0.5* 2.0 ** np.random.uniform(0,5.0)
    else:
        mu = PARAMS[0]
        sigma = mu = PARAMS[1]
        _, counts = generate_lognormal_counts(mu_nm=mu,
        sigma_nm=sigma, N_total=np.random.uniform(300, 800), CONST = CONST)
        true_counts += counts
        C_IMP = PARAMS[2]
        GAMMA0 = PARAMS[3]

    true_counts /= np.max(true_counts)  # Normalize

    omega_axis_prim = np.linspace(1280, 1360, npts)
    CONST_prim = copy.deepcopy(CONST)
    CONST_prim["OMEGA_AXIS"] = omega_axis_prim
    known_spectrum_prim = total_raman_spectrum(true_counts, C_IMP, GAMMA0, CONST=CONST_prim)

    noise_1 = np.random.uniform(-noiselevel, noiselevel, size=omega_axis_prim.shape)
    noise_2 = np.random.uniform(-noiselevel, noiselevel, size=omega_axis_prim.shape)

    omega_axis_prim = omega_axis_prim + noise_1
    known_spectrum_prim = known_spectrum_prim + noise_1
    known_spectrum_prim = known_spectrum_prim * 2.0**np.random.uniform(-4,8)
    # Merge into Nx2 array
    merged_array = np.column_stack((omega_axis_prim, known_spectrum_prim))  # Alternatively np.vstack((X_noisy, Y_noisy)).T

    # Save as CSV
    #np.savetxt("noisy_data.csv", merged_array, delimiter=",", header="X,Y", comments="", fmt="%.6f")
    np.savetxt("noisy_data.csv", merged_array, delimiter=",", comments="", fmt="%.6f")


    # Plot spectrum match
    plt.figure(figsize=(7, 4))
    plt.plot(merged_array[:,0], merged_array[:,1], label='Known spectrum', linewidth=2)
    plt.xlabel('ω (cm⁻¹)')
    plt.ylabel('I(ω)')
    plt.title('Spectrum reconstruction')
    plt.legend()
    plt.show()




    return true_counts, C_IMP, GAMMA0


#==========================================================================================



#==========================================================================================

#  T R A I N

#==========================================================================================



#==========================================================================================




def proportional_one_hot(value, length):
    # Clamp the value to the valid range [0, length-1]
    value = np.clip(value, 0, length - 1)
    
    # Find the floor and ceiling indices
    lower_idx = int(np.floor(value))
    upper_idx = int(np.ceil(value))
    
    # If value is exactly an integer, use standard one-hot
    if lower_idx == upper_idx:
        arr = np.zeros(length)
        arr[lower_idx] = 1.0
        return arr
    
    # Calculate weights based on proximity
    lower_weight = upper_idx - value  # closer to upper = less weight for lower
    upper_weight = value - lower_idx  # closer to lower = less weight for upper
    
    arr = np.zeros(length)
    arr[lower_idx] = lower_weight
    arr[upper_idx] = upper_weight
    
    return arr



def one_hot(index, length, proportional = True):
    if proportional:
        return proportional_one_hot(index, length)
    else:
        #arr = np.zeros(length)
        #arr[index] = 1.0
        #return arr
        arr = np.zeros(length)
        arr[np.rint(index).astype(int)] = 1.0
        return arr
        

def generate_sample_train(CONST, n_distributions=1, noiselevel = 0.01, proportional = True):

    # Generate composite size distribution
    total_counts = np.zeros(CONST["N_BINS"])

    if n_distributions == 2 :
         n_distributions_prim = np.random.choice(2, 1, p=[0.75,0.25])[0]+1
    elif n_distributions == 3:
        n_distributions_prim = np.random.choice(3, 1, p=[0.7,0.2,0.1])[0]+1
    else:
        n_distributions_prim = 1

    for _ in range(n_distributions_prim):
        #mu = 1.5*2.0**np.random.uniform(0.0,2.2)
        mu = np.random.uniform(1.5, 7.0)
        if np.random.uniform()>=0.32:
            #sigma = np.random.uniform(0.05, 0.5)
            sigma = np.random.uniform(0.15*(1.5/mu)**0.33, 0.4*(1.5/mu)**0.5)
            _, counts = generate_lognormal_counts(mu_nm=mu, sigma_nm=sigma, N_total=1000, CONST=CONST)
        else:
            sigma = np.random.uniform(0.35 * (mu / 1.5)**0.5, 1.3)
            _, counts = generate_normal_counts(mu_nm=mu, sigma_nm=sigma, N_total=1000, CONST=CONST)  #added Gaussian

        pow_hist = np.random.uniform(0.9,1.2)  #added deformation of distribution!
        counts2 = counts**pow_hist     #added deformation of distribution!
        counts2 = counts2 / np.sum(counts2) * np.sum(counts)     #added deformation of distribution!
        total_counts += counts2
    

    total_counts = total_counts/np.max(total_counts)
    # Choose discrete c_imp and Gamma0
    ###i_c = np.random.randint(0, CONST["N_Cimp"])
    ###i_g = np.random.randint(0, CONST["N_Gamma0"])
    i_c = np.random.uniform(0, CONST["N_Cimp"]-1)
    i_g = np.random.uniform(0, CONST["N_Gamma0"]-1)
    

    c_imp = c_imp_vs_index(i_c)
    Gamma0 = Gamma0_vs_index(i_g)

    # Generate spectrum
    spectrum = total_raman_spectrum(total_counts, c_imp, Gamma0, CONST=CONST)
    spectrum = spectrum + np.random.uniform() * np.random.uniform(-noiselevel,noiselevel,size=spectrum.shape)
    # Construct label vector Y
    total_counts_normalized = total_counts * CONST["NORMALIZE"]
    total_counts_normalized = total_counts_normalized/np.max(total_counts_normalized)
    Y = np.concatenate([
        total_counts_normalized,         # size distribution counts
        one_hot(i_c, CONST["N_Cimp"], proportional = proportional),     # one-hot c_imp
        one_hot(i_g, CONST["N_Gamma0"], proportional = proportional)       # one-hot Gamma0
    ])

    return spectrum, Y













#=====================================================================

#      C   S   V




#			  CCCCC
#			C     C
#			C
#			C
#			C
#			C     C
#			  CCCCC



#=====================================================================


import csv

def check_csv_format(filename, CONST = CONST):
    try:
        data = np.loadtxt(filename, delimiter=",")
    except Exception as e:
        return False, f"Failed to read CSV as numeric data: {e}", None

    # Check shape
    if data.ndim != 2 or data.shape[1] != 2:
        return False, "CSV must have exactly two columns", None

    X, Y = data[:, 0], data[:, 1]

    # Check positivity
    #if not np.all(X > 0) or not np.all(Y > 0):
    if not np.all(X > 0):
        return False, "All X values must be positive", None

    # Check monotonic increase
    if not np.all(np.diff(X) > 0):
        return False, "X must be strictly increasing", None

    # Find s and f indices
    s_candidates = np.where((X[:-1] < 1300) & (X[1:] > 1300))[0]
    f_candidates = np.where((X[:-1] < 1340) & (X[1:] > 1340))[0]
    if len(s_candidates) == 0 or len(f_candidates) == 0:
        return False, "Cannot find interpolation bounds for 1300 and 1350", None

    s = s_candidates[0]
    f = f_candidates[0]

    if f <= s:
        return False, "f must be greater than s", None

    sub_X = X[s:f+2]
    sub_Y = Y[s:f+2]

    if len(sub_X) < 40:
        return False, f"Subarray has only {len(sub_X)} points, at least 40 required", None

    # Check uniformity
    deltas = np.diff(sub_X)
    mean_delta = np.mean(deltas)
    if not np.all((deltas > 0.01 * mean_delta) & (deltas < 10 * mean_delta)):
        return False, "Non-uniform spacing in X subarray", None

    # Target X values for interpolation
    #interp_X = np.arange(1300, 1340.1, 0.5)
    interp_X = CONST["OMEGA_AXIS"]
    # Homemade linear interpolation
    interp_Y = np.empty_like(interp_X)
    j = 0
    for i, xi in enumerate(interp_X):
        # Move j to correct interval
        while j < len(sub_X) - 2 and sub_X[j+1] < xi:
            j += 1
        x0, x1 = sub_X[j], sub_X[j+1]
        y0, y1 = sub_Y[j], sub_Y[j+1]
        # Linear interpolation formula
        t = (xi - x0) / (x1 - x0)
        interp_Y[i] = (1 - t) * y0 + t * y1
    interp_Y = interp_Y / np.max(interp_Y)
    interpolated_data = np.column_stack((interp_X, interp_Y))
    return True, "Success", interpolated_data




#=====================================================================

#   F I G U R E S

#   for calculation of spectrum


#			FFFFFFF
#			F
#			F
#			FFFF
#			F
#			F
#			F




#=====================================================================









def create_output_zip_derive(zip_path, input1, CONST = CONST, Cleanup = True):
    # Prepare paths
    npts = 100;
    omega_axis_prim = np.linspace(1280, 1360, npts)
    CONST_prim = copy.deepcopy(CONST)
    CONST_prim["OMEGA_AXIS"] = omega_axis_prim

    P1, P2, P3, P4 = input1;

    _, derived_counts =  generate_lognormal_counts(mu_nm=P1, sigma_nm=P2, N_total=1000, CONST = CONST_prim)
    derived_counts = derived_counts/np.max(derived_counts)
    derived_spectrum = total_raman_spectrum(derived_counts, P3, P4, n_max=5, CONST = CONST_prim)

    ###base_dir = os.getcwd()
    base_dir = os.path.dirname(zip_path)


    spectra_csv_path = os.path.join(base_dir, "spectrum.csv")
    hist_csv_path = os.path.join(base_dir, "histogram.csv")

    spectra_plot_path = os.path.join(base_dir, "spectrum.png")
    hist_plot_path = os.path.join(base_dir, "histogram.png")

    # Save processed CSV
    exp_ram = np.column_stack((CONST_prim["OMEGA_AXIS"],derived_spectrum))
    np.savetxt(spectra_csv_path, exp_ram, delimiter=",", fmt="%.6f")
    exp_hist = np.column_stack((CONST_prim["BINS"],derived_counts))
    np.savetxt(hist_csv_path, exp_hist, delimiter=",", fmt="%.6f")

    #==========================================================================


    # Plot spectrum match
    plt.figure(figsize=(7, 4))
    plt.plot(exp_ram[:,0], exp_ram[:,1], label='Derived spectrum', linewidth=2)
    plt.xlabel('ω (cm⁻¹)')
    plt.ylabel('I(ω)')
    plt.title("Custom Raman spectrum")
    plt.legend()
    plt.savefig(spectra_plot_path)
    plt.close()

    # Plot distribution match
    plt.figure(figsize=(7, 4))
    #plt.bar(sizes_nm - 0.1, true_counts, width=0.2, label='True distribution')
    #plt.bar(sizes_nm + 0.1, reconstructed_distribution, width=0.2, label='Reconstructed distribution')
    plt.bar(exp_hist[:,0], exp_hist[:,1], width=0.29, label='Reconstructed distribution')
    plt.xlabel('Particle size (nm)')
    plt.ylabel('Normalized counts')
    plt.title(f"Size distribution: $C_{{\\mathrm{{imp}}}} = {P3:.2f}$%  $\Gamma_0 = {P4:.2f}~\\mathrm{{cm}}^{{-1}}$")
    plt.legend()
    plt.savefig(hist_plot_path)
    plt.close()

    #==========================================================================

    # Create ZIP archive
    with zipfile.ZipFile(zip_path, 'w') as zipf:
        zipf.write(spectra_csv_path, arcname="spectrum.csv")
        zipf.write(hist_csv_path, arcname="histogram.csv")

        zipf.write(spectra_plot_path, arcname="spectrum.png")
        zipf.write(hist_plot_path, arcname="histogram.png")

    # Cleanup
    if Cleanup:
        os.remove(spectra_csv_path)
        os.remove(hist_csv_path)
        os.remove(spectra_plot_path)
        os.remove(hist_plot_path)






#====================================================================================================
#      N E U R A L


#			M      M
#			MM     M
#			M M    M
#			M  M   M
#			M    M M
#			M     MM
#			M      M


#====================================================================================================





def load_numpy_model(filename="model_weights.npz", activation='tanh'):
    data = np.load(filename)
    weights = [data[key] for key in data]

    layers = []
    for i in range(0, len(weights), 2):
        W = weights[i]
        b = weights[i+1]
        layers.append((W, b))

    if activation == 'relu':
        act_fn = lambda x: np.maximum(0, x)
    elif activation == 'tanh':
        act_fn = np.tanh
    else:
        raise ValueError("Only 'relu' and 'tanh' are supported.")

    def model(x):
        out = x
        for i, (W, b) in enumerate(layers):
            out = out @ W.T + b
            if i < len(layers) - 1:
                out = act_fn(out)
        return out

    return model

def weighted_round_index(arr, pow1 = 1.0):
    indices = np.arange(arr.shape[0])
    if np.sum(arr) == 0:  # avoid divide by zero
        return int(np.round(np.mean(indices)))
    mean_index = np.sum((indices**pow1) * arr) / np.sum(arr)
    #return int(np.round(mean_index))
    return mean_index










def NNapproach(known_spectrum, CONST = CONST, PARAMS = False):
    spectrum_to_treat = known_spectrum[:,1]
    spectrum_to_treat = spectrum_to_treat + np.random.uniform(low=-0.005*0, high=0.005*3, size=spectrum_to_treat.shape)
    spectrum_to_treat = spectrum_to_treat/np.max(spectrum_to_treat)

    base_dir_NN = os.path.dirname(os.path.abspath(__file__))
    if PARAMS == False:
        wts_NN = os.path.join(base_dir_NN, "model_weights_here.npz")
    else:
        wts_NN = os.path.join(base_dir_NN, PARAMS[0])

    numpy_model = load_numpy_model(wts_NN, activation='tanh')
    Y_pred = numpy_model(spectrum_to_treat)
    Y_pred= np.clip(Y_pred,a_min=0.0,a_max=2.0)

    bins_return = Y_pred[:CONST["N_BINS"]]
    bins_return = np.clip(bins_return,a_min=0.0, a_max=10000)
    bins_return = bins_return / CONST["NORMALIZE"]
    bins_return = bins_return/np.max(bins_return)
    #------------------------
    

    cimp_arr_tmp = Y_pred[CONST["N_BINS"]:CONST["N_BINS"] + CONST["N_Cimp"]]
    gamma_arr_tmp = Y_pred[-CONST["N_Gamma0"]:]
    i_c_pred = np.argmax(cimp_arr_tmp)
    i_g_pred = np.argmax(gamma_arr_tmp)
    c_imp_pred = c_imp_vs_index(i_c_pred)
    Gamma0_pred = Gamma0_vs_index(i_g_pred)
    #------------------------
    if not (PARAMS==False):
        if PARAMS[1] == 1:
            current_msd = 100000000.0
            argsort_cimp_arr = np.argsort(cimp_arr_tmp)[-5:][::-1]
            for i_c_tmp in argsort_cimp_arr:
                trial_spectrum = total_raman_spectrum(
                    distribution_counts=bins_return,
                    c_imp=c_imp_vs_index(i_c_tmp),
                    Gamma0=Gamma0_pred,
                    CONST=CONST
                )
                trial_msd = mean_squared_deviation(trial_spectrum, known_spectrum[:,1])
                delta_msd = trial_msd - current_msd
                if delta_msd < 0:
                    i_c_pred = i_c_tmp
                    current_msd = trial_msd
    c_imp_pred = c_imp_vs_index(i_c_pred)
    #-------------------------


    spectrum_pred = total_raman_spectrum(
        distribution_counts=bins_return,
        c_imp=c_imp_pred,
        Gamma0=Gamma0_pred,
        CONST=CONST
    )


    result_tmp = known_spectrum, np.column_stack((CONST["BINS"],bins_return)), np.column_stack((CONST["OMEGA_AXIS"],spectrum_pred)), c_imp_pred, Gamma0_pred, Y_pred
    return result_tmp






