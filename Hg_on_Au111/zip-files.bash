#!/bin/bash

# Define the target zip archive name
ZIP_NAME="Hg_on_Au111_modular.zip"

# Define the list of files required for the modular workflow
FILES=(
    "calculator.py"
    "cli.py"
    "config.py"
    "constants.py"
    "energies.py"
    "optimization.py"
    "phonons.py"
    "run_hg_au111.py"
    "surface.py"
    "thermodynamics.py"
    "visualize_hg_vibrations.py"
    "workflow.py"
    "requirements.txt"
    "hg_au111_config.ini"
    "Hg_on_Au111_Theory_and_Code_Description.rst"
)

echo "Verifying file existence before archiving..."
MISSING_FILE=0

for file in "${FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Error: Required file '$file' is missing."
        MISSING_FILE=1
    fi
done

if [ $MISSING_FILE -eq 1 ]; then
    echo "Archive creation aborted due to missing files."
    exit 1
fi

echo "All files verified. Archiving into $ZIP_NAME..."
zip "$ZIP_NAME" "${FILES[@]}"

if [ $? -eq 0 ]; then
    echo "Success! $ZIP_NAME has been created successfully."
else
    echo "An error occurred during the zip creation process."
    exit 1
fi

