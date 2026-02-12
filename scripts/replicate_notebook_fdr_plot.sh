#!/usr/bin/env bash
# Replicate the FDR validation plot from the original aa.fdr_github.ipynb notebook
#
# This script runs the FDR validation using the ABRF dataset that can be downloaded
# with: python download_data.py
#
# The notebook uses de novo results from different "samples" as stand-in decoys
# to demonstrate how decoy-based FDR estimation compares to true (database) FDR.
# Samples 3-7 are used as "decoys" against Sample 10 as the target.
#
# Usage:
#   ./scripts/replicate_notebook_fdr_plot.sh
#
# Output:
#   fig/notebook_fdr_validation.png

set -e

# Path configuration - adjust if your data is in a different location
DATA_DIR="data"
OUTPUT_DIR="fig"

# Input files from the ABRF dataset
SPECTRUM_FILE="${DATA_DIR}/2017-12-4_ABRF_200_DDA1.mgf"
DB_FILE="${DATA_DIR}/pd_merged.csv.db.psms.csv"
TARGET_FILE="${DATA_DIR}/PEAKS/Sample 10.denovo.csv"

# Decoy files (using samples 3-7 as in the notebook)
# The notebook labels these as "decoy 70%", "decoy 60%", etc.
DECOY_FILES=(
    "${DATA_DIR}/PEAKS/Sample 3.denovo.csv"
    "${DATA_DIR}/PEAKS/Sample 4.denovo.csv"
    "${DATA_DIR}/PEAKS/Sample 5.denovo.csv"
    "${DATA_DIR}/PEAKS/Sample 6.denovo.csv"
    "${DATA_DIR}/PEAKS/Sample 7.denovo.csv"
)

# Labels matching the notebook (10 - sample_number)0%
LABELS=(
    "decoy 70%"
    "decoy 60%"
    "decoy 50%"
    "decoy 40%"
    "decoy 30%"
)

OUTPUT_FILE="${OUTPUT_DIR}/notebook_fdr_validation.png"

# Check that required files exist
echo "Checking input files..."
for file in "$SPECTRUM_FILE" "$DB_FILE" "$TARGET_FILE" "${DECOY_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file not found: $file"
        echo ""
        echo "Please download the data first:"
        echo "  python download_data.py"
        exit 1
    fi
done

# Create output directory if needed
mkdir -p "$OUTPUT_DIR"

echo ""
echo "Running FDR validation..."
echo "  Target: $TARGET_FILE"
echo "  Decoys: ${#DECOY_FILES[@]} files"
echo "  Database: $DB_FILE"
echo "  Spectra: $SPECTRUM_FILE"
echo "  Output: $OUTPUT_FILE"
echo ""

# Run the FDR validation command
# Note: Using --fdr-max 0.05 to match the notebook's xlim/ylim settings
novoboard fdr \
    --target-file "$TARGET_FILE" \
    --decoy-files "${DECOY_FILES[@]}" \
    --db-file "$DB_FILE" \
    --spectrum-file "$SPECTRUM_FILE" \
    --output-file "$OUTPUT_FILE" \
    --score-column "ALC (%)" \
    --aa-score-column "local confidence (%)" \
    --ion-threshold 0.90 \
    --labels "${LABELS[@]}" \
    --fdr-max 0.05 \
    --dpi 150

echo ""
echo "Done! Plot saved to: $OUTPUT_FILE"
