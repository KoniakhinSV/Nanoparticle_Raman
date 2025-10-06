import raman_routines
import os

# --- initialize ---
os.environ["QT_QPA_PLATFORM"] = "offscreen"

CONST = raman_routines.define_constans()
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
original_path = os.path.join(BASE_DIR, "input.csv")
zip_path = os.path.join(BASE_DIR, "result_Raman_Fit_NN.zip")

#process
is_valid, message, known_spectrum = raman_routines.check_csv_format(original_path, CONST=CONST)
if not is_valid:
    print(f"Invalid CSV: {message}")
else:
    processing_function = raman_routines.NNapproach
    params1 = [raman_routines.find_load_npz_file_alt(), 1]
    result_main = processing_function(known_spectrum=known_spectrum, CONST=CONST,PARAMS=params1)
    raman_routines.create_output_zip(
        original_path,
        zip_path,
        result_main,
        mode = "NN"
    )
