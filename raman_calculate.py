import raman_routines
from csv_param import validate_parameters
import os
import csv

#os.environ["QT_QPA_PLATFORM"] = "offscreen"

CONST = raman_routines.define_constans()
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
zip_path = os.path.join(BASE_DIR, "result_Raman_Calculated.zip")

params, error = validate_parameters('parameters.csv')
if error:
    print(f"Validation error: {error}")
else:
    print(f"Read parameters are: {params}")
    param_list = [params["D"],params["dD"],params["Cimp"],params["Gamma0"]]
    raman_routines.create_output_zip_derive(zip_path, param_list, CONST = CONST)