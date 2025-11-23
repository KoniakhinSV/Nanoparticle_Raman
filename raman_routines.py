import numpy as np
import matplotlib
matplotlib.use("Agg")    # uncomment for standalone
import matplotlib.pyplot as plt

import os
import zipfile
import copy
import csv


#from site 31 OCT 2025 ===================================================================

def find_load_npz_file_alt(endswith_argument = '_load.npz'):

    #cwd = os.getcwd()
    cwd = os.path.dirname(os.path.abspath(__file__))

    for filename in os.listdir(cwd):
        if filename.endswith(endswith_argument):
            return os.path.join(os.path.dirname(os.path.abspath(__file__)),filename)


#===================================================================


def estimate_normal_params(x, y):
    """
    Estimate mean and standard deviation from histogram data
    
    Parameters:
    x: bin centers or bin edges
    y: bin counts or frequencies
    
    Returns:
    mean, std: estimated parameters
    """
    # If x represents bin edges, convert to bin centers
    if len(x) == len(y) + 1:
        x_centers = (x[:-1] + x[1:]) / 2
    else:
        x_centers = x
    
    # Normalize y to get probabilities
    probabilities = y / np.sum(y)
    
    # Calculate mean
    mean = np.average(x_centers, weights=probabilities)
    
    # Calculate standard deviation
    variance = np.average((x_centers - mean)**2, weights=probabilities)
    std = np.sqrt(variance)
    
    return mean, std



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


    # --- Step 1: Initialize current distribution from natural ND parameters / take initial params ---
    C_IMP = 0.0
    GAMMA0 = 15.0
    #mu0 = np.random.uniform(2.0, 7.0)
    #sigma0 = np.random.uniform(0.2, 0.6)
    mu0 = 4.0
    sigma0 = 0.3
    Adjust_C_IMP = 1.0
    Adjust_GAMMA0 = 1.0
    TEMP0 = 0.000001
    if PARAMS and isinstance(PARAMS, (list, tuple)):
        if len(PARAMS) == 6:
            mu0, sigma0, C_IMP, GAMMA0, Adjust_C_IMP, Adjust_GAMMA0 = PARAMS
        elif len(PARAMS) == 7:
            mu0, sigma0, C_IMP, GAMMA0, Adjust_C_IMP, Adjust_GAMMA0, TEMP0 = PARAMS


    _, current_counts = generate_lognormal_counts(mu_nm=mu0, sigma_nm=sigma0, N_total=1000, CONST = CONST)
    current_counts /= np.max(current_counts)
    current_spectrum = total_raman_spectrum(current_counts, C_IMP, GAMMA0, CONST=CONST)
    current_msd = mean_squared_deviation(current_spectrum, known_spectrum[:,1])

    # --- Step 2: Iterative Metropolis-like refinement ---

    history_msd = [current_msd]
    history_counts = [current_counts.copy()]
    MAX_ITER = 1000 * 3

    for step in range(MAX_ITER):
        # Generate small correction
        mu_c = np.random.uniform(2.0, 7.0)
        sigma_c = np.random.uniform(0.05, 0.3*6/mu_c)
        _, correction = generate_lognormal_counts(mu_c, sigma_c, N_total=1000, CONST=CONST)
        correction /= np.max(correction)

        # Add correction
        trial_counts = current_counts + 0.05 * correction * np.random.rand()
        #trial_counts = current_counts + 0.05 * correction * (2*np.random.rand() - 1)
        trial_counts = np.clip(trial_counts,
                                a_min = 0.0, a_max = 40.0);
        trial_counts /= np.max(trial_counts)  # keep normalized

        C_IMP_trial = np.clip(C_IMP + Adjust_C_IMP * np.random.uniform(-0.01,0.01),
                                a_min = 0.0, a_max = 3.0);
        GAMMA0_trial = np.clip(GAMMA0 + Adjust_GAMMA0* np.random.uniform(-0.05,0.05),
                                a_min = 0.5, a_max = 40.0);

        trial_spectrum = total_raman_spectrum(trial_counts, C_IMP_trial, GAMMA0_trial,
        CONST=CONST)
        trial_msd = mean_squared_deviation(trial_spectrum, known_spectrum[:,1])
        delta_msd = trial_msd - current_msd
        rr = np.random.uniform()
        if np.exp(-delta_msd/TEMP0) > rr: ##uncomment for full Metropolis
        #if delta_msd < 0:   # comment for full Metropolis
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

    #mean_D = np.sum(reconstructed_distribution[:,0] * reconstructed_distribution[:,1])/np.sum(reconstructed_distribution[:,1])
    #mean_SCA = np.sum(reconstructed_distribution[:,0]**4 * reconstructed_distribution[:,1])/np.sum(reconstructed_distribution[:,0]**3 * reconstructed_distribution[:,1])

    mean_N, sigma_N = estimate_normal_params(reconstructed_distribution[:,0], reconstructed_distribution[:,1])
    mean_SCA, sigma_SCA = estimate_normal_params(reconstructed_distribution[:,0], reconstructed_distribution[:,1]*(reconstructed_distribution[:,0]**3))

    base_dir = os.path.dirname(original_csv_path)
    original_csv_name = os.path.basename(original_csv_path)

    spectra_csv_path = os.path.join(base_dir, mode+"_spectra.csv")
    hist_csv_path = os.path.join(base_dir, mode+"_histogram.csv")

    spectra_plot_path = os.path.join(base_dir, mode+"_spectra.png")
    hist_plot_path = os.path.join(base_dir, mode+"_histogram.png")

    # Save processed CSV

    if known_spectrum.shape[1] == 3:
        exp_ram = np.column_stack((known_spectrum[:,0], known_spectrum[:,1], known_spectrum[:,2], reconstructed_spectrum[:,1]))
        np.savetxt(spectra_csv_path, exp_ram, delimiter=",", fmt="%.6f", header='Energy, Treated, Original, Reconstructed')
    else:
        exp_ram = np.column_stack((known_spectrum[:,0], known_spectrum[:,1],reconstructed_spectrum[:,1]))
        np.savetxt(spectra_csv_path, exp_ram, delimiter=",", fmt="%.6f", header='Energy, Treated, Reconstructed')

    #exp_hist = np.column_stack((reconstructed_distribution))
    exp_hist = np.column_stack((reconstructed_distribution[:,0],reconstructed_distribution[:,1],volume_distribution))
    np.savetxt(hist_csv_path, exp_hist, delimiter=",", fmt="%.6f", header='Size_nm, By_N, By_Vol')

    #==========================================================================


    # Plot spectrum match
    plt.figure(figsize=(10, 6))
    plt.plot(known_spectrum[:,0], known_spectrum[:,1], label='Treated spectrum', linewidth=2, color = 'black')
    plt.plot(known_spectrum[:,0], reconstructed_spectrum[:,1], label='Reconstructed spectrum', linestyle='--',color = 'blue')

    if known_spectrum.shape[1] == 3:
        plt.plot(known_spectrum[:,0], known_spectrum[:,2], label='With background (original)', color='green', linestyle='dotted')

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
    plt.text(0.0,-0.13,f"$D_N = ${mean_N:.2f}nm  $D_V = ${mean_SCA:.2f}nm  $\\Delta D_V = ${sigma_SCA:.2f}nm")

    plt.legend()
    plt.savefig(hist_plot_path)
    plt.close()

    #==========================================================================

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
    plt.bar(exp_hist[:,0], exp_hist[:,1], width=0.29, label='Derived distribution')
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
    spectrum_to_treat = spectrum_to_treat + np.random.uniform(low=-0.005*0, high=0.004, size=spectrum_to_treat.shape)
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




#=============================================================================================================================================










#==============================================================================================================================================

#      C   S   V




#			  CCCCC
#			C     C
#			C
#			C
#			C
#			C     C
#			  CCCCC



#=====================================================================

def subtract_background(known_spectrum_arg, filename_BACK):
    known_spectrum = known_spectrum_arg.copy()
    numpy_model_BACK = load_numpy_model(filename =  filename_BACK , activation='tanh')
    known_spectrum_BACK = numpy_model_BACK(known_spectrum[:,1])
    known_spectrum[:,1] = known_spectrum_BACK[:]
    sigma = 2 * 0.67
    radius = int(3 * sigma)
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (x / sigma) ** 2)
    kernel = kernel / kernel.sum()
    pad_width = radius
    padded_data = np.pad(known_spectrum[:,1], pad_width, mode='edge')
    known_spectrum[:,1] = np.convolve(padded_data, kernel, mode='valid')
    #known_spectrum[:,1] = np.convolve(known_spectrum[:,1], kernel, mode='same')
    known_spectrum[:,1] = known_spectrum[:,1] / np.max(known_spectrum[:,1])
    return np.column_stack((known_spectrum[:,0], known_spectrum[:,1],known_spectrum_arg[:,1]))

#=====================================================================







def check_csv_format(filename, CONST=None, SUBTRACT_BACKGROUND=False):
    if CONST is None:
        CONST = {"OMEGA_AXIS": np.arange(1300, 1340.1, 0.5)}

    try:
        # First, try to detect the format by reading the first few lines
        with open(filename, 'r', encoding='utf-8') as file:
            lines = [file.readline().strip() for _ in range(5)]

        # Detect delimiter
        delimiter = detect_delimiter(lines[0])

        # Detect decimal separator and check for header
        has_header, decimal_separator = detect_header_and_decimal(lines, delimiter)

        # Load data, skipping header if present
        skiprows = 1 if has_header else 0

        if decimal_separator == ',':
            # For comma decimal, we need to convert commas to dots
            data = load_csv_with_comma_decimal(filename, delimiter, skiprows)
        else:
            # Standard load with detected delimiter
            data = np.loadtxt(filename, delimiter=delimiter, skiprows=skiprows)

    except Exception as e:
        return False, f"Failed to read CSV file: {e}", None

    # Check shape
    if data.ndim != 2 or data.shape[1] != 2:
        return False, f"CSV must have exactly two columns, found {data.shape[1] if data.ndim == 2 else 'invalid'}", None

    X, Y = data[:, 0], data[:, 1]

    # Check positivity
    if not np.all(X > 0):
        return False, "All X values must be positive", None

    # Check monotonic increase
    if not np.all(np.diff(X) > 0):
        return False, "X must be strictly increasing", None

    # Find s and f indices
    s_candidates = np.where((X[:-1] <= 1300) & (X[1:] > 1300))[0]
    f_candidates = np.where((X[:-1] < 1340) & (X[1:] >= 1340))[0]
    if len(s_candidates) == 0 or len(f_candidates) == 0:
        return False, "Cannot find interpolation bounds for 1300 and 1340", None

    s = s_candidates[0]
    f = f_candidates[0]

    if f <= s:
        return False, "f must be greater than s", None

    sub_X = X[s:f+2]
    sub_Y = Y[s:f+2]

    if len(sub_X) < 20:
        return False, f"Subarray has only {len(sub_X)} points, at least 20 required", None

    # Check uniformity
    deltas = np.diff(sub_X)
    mean_delta = np.mean(deltas)
    if not np.all((deltas > 0.01 * mean_delta) & (deltas < 10 * mean_delta)):
        return False, "Non-uniform spacing in X subarray", None

    # Target X values for interpolation
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

    if SUBTRACT_BACKGROUND:
        filename_BACK = find_load_npz_file_alt(endswith_argument="loadB.npz")
        if filename_BACK is not None:
            interpolated_data = subtract_background(interpolated_data, filename_BACK)

    return True, "Success", interpolated_data


def detect_delimiter(first_line):
    """Detect the delimiter used in the CSV file"""
    # Common delimiters to check
    delimiters = [',', '\t', ';', ' ', '|']

    # Count occurrences of each delimiter
    counts = {delim: first_line.count(delim) for delim in delimiters}

    # Return the delimiter with the highest count
    best_delimiter = max(counts, key=counts.get)

    # If no delimiter found or space is the best but there are few spaces,
    # default to comma
    if counts[best_delimiter] == 0 or (best_delimiter == ' ' and counts[' '] < 2):
        return ','

    return best_delimiter


def detect_header_and_decimal(lines, delimiter):
    """Detect if the file has a header and what decimal separator is used"""
    has_header = False
    decimal_separator = '.'

    if not lines or not lines[0]:
        return has_header, decimal_separator

    # Check if first line might be a header (contains non-numeric values)
    first_line_parts = lines[0].split(delimiter)
    try:
        # Try to convert first two elements to float
        float(first_line_parts[0].strip())
        float(first_line_parts[1].strip())
        # If successful, likely no header
        has_header = False
    except (ValueError, IndexError):
        # If conversion fails, likely has header
        has_header = True

    # Detect decimal separator by checking numeric lines
    check_lines = lines[1:] if has_header else lines

    for line in check_lines:
        if not line:
            continue

        parts = line.split(delimiter)
        if len(parts) >= 2:
            # Check for comma as decimal separator
            for part in parts[:2]:  # Check first two columns
                part = part.strip()
                if ',' in part and '.' not in part:
                    # Count commas - if exactly one and not at the end, likely decimal
                    if part.count(',') == 1 and not part.endswith(','):
                        try:
                            # Try replacing comma with dot and converting
                            float(part.replace(',', '.'))
                            decimal_separator = ','
                            return has_header, decimal_separator
                        except ValueError:
                            pass
                elif '.' in part:
                    decimal_separator = '.'

    return has_header, decimal_separator


def load_csv_with_comma_decimal(filename, delimiter, skiprows):
    """Load CSV file with comma as decimal separator"""
    def convert_comma_decimal(x):
        if isinstance(x, bytes):
            x = x.decode('utf-8')
        return float(x.replace(',', '.'))

    converters = {0: convert_comma_decimal, 1: convert_comma_decimal}
    return np.loadtxt(filename, delimiter=delimiter, skiprows=skiprows, converters=converters)

