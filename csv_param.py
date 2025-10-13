import raman_routines
import os
import csv


#=====================================================================================
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
                    'background': float(row.get('background')),
                    'D': float(row.get('D')),
                    'dD': float(row.get('dD')),
                    'Cimp': float(row.get('Cimp')),
                    'Gamma0': float(row.get('Gamma0')),
                    "adjust_Cimp": float(row.get('adjust_Cimp')),
                    "adjust_Gamma0": float(row.get('adjust_Gamma0'))
                }
            except (ValueError, TypeError) as e:
                return None, f"Invalid parameter values in CSV: {str(e)}"

            
            if not (params['background'] == 0 or params['background'] == 1):
                return None, f"background parameter 1.0 or 0.0, got {params['background']}"
            if not (params['adjust_Cimp'] == 0 or params['adjust_Cimp'] == 1):
                return None, f"adjust_Cimp parameter 1.0 or 0.0, got {params['adjust_Cimp']}"
            if not (params['adjust_Gamma0'] == 0 or params['adjust_Gamma0'] == 1):
                return None, f"backgradjust_Gamma0oung parameter 1.0 or 0.0, got {params['adjust_Gamma0']}"


            # =================================================
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
            

            return params, None

    except FileNotFoundError:
        return None, f"File not found: {csv_file_path}"
    except Exception as e:
        return None, f"Error reading CSV file: {str(e)}"
#=====================================================================================