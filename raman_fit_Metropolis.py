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
params, error = raman_routines.validate_parameters('parameters.csv')
if error:
    print(f"Parameters validation error: {error}")
    params1 = False

else:
    print(f"Read initial parameters are: {params}")
    params1 = [params["D"],params["dD"],params["Cimp"],params["Gamma0"],params["adjust_Cimp"],params["adjust_Gamma0"]]
#----------------------------------------------------------------------------------

    #process
    if params["background"] == 1.0:
        SUBTRACT_BACKGROUND = True
    else:
        SUBTRACT_BACKGROUND = False

    is_valid, message, known_spectrum = raman_routines.check_csv_format(original_path, CONST=CONST, SUBTRACT_BACKGROUND = SUBTRACT_BACKGROUND)
    if not is_valid:
        print(f"Invalid CSV with spectrum: {message}")
    else:
        processing_function = raman_routines.metropolis
        result_main = processing_function(known_spectrum=known_spectrum, CONST = CONST, PARAMS = params1)
        raman_routines.create_output_zip(
            original_path,
            zip_path,
            result_main,
            mode = "M"
        )
