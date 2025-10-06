import raman_routines
import os
import csv

# --- initialize ---
os.environ["QT_QPA_PLATFORM"] = "offscreen"

CONST = raman_routines.define_constans()
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
original_path = os.path.join(BASE_DIR, "input.csv")
zip_path = os.path.join(BASE_DIR, "result_Raman_Fit_Metropolis.zip")


#----------------------------------------------------------------------------------

def validate_parameters(csv_file_path):
    """
    Reads and validates parameters from a CSV file according to specified constraints.

    Args:
        csv_file_path (str): Path to the CSV file containing parameters

    Returns:
        dict: Dictionary with parameter values if valid, None if invalid
        str: Error message if validation fails, None otherwise
    """
    try:
        with open(csv_file_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            row = next(reader)  # Read the first row

            # Extract parameters
            try:
                params = {
                    'D': float(row.get('D')),
                    'dD': float(row.get('dD')),
                    'Cimp': float(row.get('Cimp')),
                    'Gamma0': float(row.get('Gamma0'))
                }
            except (ValueError, TypeError) as e:
                return None, f"Invalid parameter values in CSV: {str(e)}"

            # Validate P1: Size D_0 (1.5-9)
            if not (1.5 <= params['D'] <= 9):
                return None, f"Diameter D must be between 1.5 and 9, got {params['D']}"

            # Validate P2: Size dispersion (0-1)
            if not (0.05 <= params['dD'] <= 1.0):
                return None, f"Size dispersion dD must be between 0.05 and 1.0, got {params['P2']}"

            # Validate P3: C_imp (0-5)
            if not (0.0 <= params['Cimp'] <= 5.0):
                return None, f"C_imp must be between 0.0 and 5.0, got {params['Cimp']}"

            # Validate P4: Gamma_0 (0-10)
            if not (0.5 <= params['Gamma0'] <= 15.0):
                return None, f"Gamma0 must be between 0.5 and 15.0, got {params['Gamma0']}"
            if 1==2:
            # Check if values are multiples of their step sizes
                if round(params['D'] * 10) % 1 != 0:
                    return None, f"D must be a multiple of 0.1, got {params['D']}"

                if round(params['dD'] * 100) % 5 != 0:
                    return None, f"dD must be a multiple of 0.05, got {params['dD']}"

                if round(params['Cimp'] * 10) % 1 != 0:
                    return None, f"Cimp must be a multiple of 0.1, got {params['Cimp']}"

                if round(params['Gamma0'] * 10) % 1 != 0:
                    return None, f"Gamma0 must be a multiple of 0.1, got {params['Gamma0']}"

            return params, None

    except FileNotFoundError:
        return None, f"File not found: {csv_file_path}"
    except Exception as e:
        return None, f"Error reading parameter CSV file: {str(e)}"


#----------------------------------------------------------------------------------
params, error = validate_parameters('parameters.csv')
if error:
    print(f"Validation error: {error}")
    params1 = False

else:
    print(f"Read initial parameters are: {params}")
    params1 = [params["D"],params["dD"],params["Cimp"],params["Gamma0"]]
#----------------------------------------------------------------------------------

#process
is_valid, message, known_spectrum = raman_routines.check_csv_format(original_path, CONST=CONST)
if not is_valid:
    print(f"Invalid CSV: {message}")
else:
    processing_function = raman_routines.metropolis
    result_main = processing_function(known_spectrum=known_spectrum, CONST = CONST, PARAMS = params1)
    raman_routines.create_output_zip(
        original_path,
        zip_path,
        result_main,
        mode = "M"
    )
