#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Path to the app (relative to the script directory)
#APP_PATH="$SCRIPT_DIR/raman_calculate.py"
#APP_PATH="$SCRIPT_DIR/raman_fit_NN.py"
APP_PATH="$SCRIPT_DIR/raman_fit_Metropolis.py"

# Function to find Python interpreter in subdirectories
find_python_interpreter() {
    local search_dir="$1"

    # Look for common Python locations in subdirectories
    local possible_paths=(
        "$search_dir"/*/bin/python
        "$search_dir"/*/bin/python3
        "$search_dir"/*/Scripts/python.exe
        "$search_dir"/*/Scripts/python3.exe
        "$search_dir"/venv/bin/python
        "$search_dir"/venv/bin/python3
        "$search_dir"/.venv/bin/python
        "$search_dir"/.venv/bin/python3
        "$search_dir"/env/bin/python
        "$search_dir"/env/bin/python3
        "$search_dir/.conda/bin/python"
        "$search_dir/.conda/bin/python3"
    )

    for path in "${possible_paths[@]}"; do
        # Expand glob pattern and check if file exists and is executable
        for expanded_path in $path; do
            if [[ -f "$expanded_path" && -x "$expanded_path" ]]; then
                echo "$expanded_path"
                return 0
            fi
        done
    done

    return 1
}

# Try to find Python interpreter in script directory
PYTHON_PATH=$(find_python_interpreter "$SCRIPT_DIR")

if [[ -n "$PYTHON_PATH" ]]; then
    echo "Found Python interpreter: $PYTHON_PATH"
    "$PYTHON_PATH" "$APP_PATH" "$@"
else
    echo "Error: No Python interpreter found in subdirectories of $SCRIPT_DIR"
    echo "Trying system Python as fallback..."

    # Fallback to system Python
    if command -v python3 &> /dev/null; then
        python3 "$APP_PATH" "$@"
    elif command -v python &> /dev/null; then
        python "$APP_PATH" "$@"
    else
        echo "Error: No Python interpreter found anywhere!"
        exit 1
    fi
fi
