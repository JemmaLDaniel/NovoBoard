#!/bin/bash
# Evaluate FDR for per-origin data using the combined predictions CSV.
# This script filters the combined predictions by experiment_name, preprocesses,
# and runs novoboard fdr for each origin.
#
# Prerequisites:
#   - Combined predictions CSV: per_origin_preds.csv
#   - Per-origin MGF files in per_origin_mgf/

set -e

# Configuration
COMBINED_PREDS="/home/j-daniel/Documents/winnow/novoboard/benchmark_data/instanovo_outputs/per_origin_preds.csv"
MGF_DIR="/home/j-daniel/Documents/winnow/novoboard/benchmark_data/instanovo_inputs/per_origin_mgf"
WORK_DIR="/home/j-daniel/Documents/winnow/novoboard/benchmark_data/per_origin_processed"
OUTPUT_DIR="/home/j-daniel/repos/NovoBoard/fig/per_origin"

# Sampling rates to evaluate (fraction of peaks KEPT)
# E.g., 0.10 = 10% kept = 90% masked, 0.90 = 90% kept = 10% masked
SAMPLING_RATES="0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90"

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$WORK_DIR"

# Check that combined predictions file exists
if [ ! -f "$COMBINED_PREDS" ]; then
    echo "Error: Combined predictions file not found: $COMBINED_PREDS"
    exit 1
fi

# Get list of origins from MGF files (excluding decoy files)
ORIGINS=$(ls "$MGF_DIR"/*.mgf 2>/dev/null | xargs -n1 basename | grep -v '_decoy_' | sed 's/.mgf$//')

if [ -z "$ORIGINS" ]; then
    echo "Error: No MGF files found in $MGF_DIR"
    exit 1
fi

echo "Found origins: $ORIGINS"
echo "Combined predictions: $COMBINED_PREDS"
echo ""

# Python script to filter CSV by experiment_name
filter_csv() {
    local input_csv="$1"
    local output_csv="$2"
    local experiment_name="$3"

    python3 << EOF
import csv
import sys

input_file = "$input_csv"
output_file = "$output_csv"
experiment_name = "$experiment_name"

with open(input_file, 'r', newline='') as infile:
    reader = csv.DictReader(infile)
    fieldnames = reader.fieldnames

    with open(output_file, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        count = 0
        for row in reader:
            if row.get('experiment_name') == experiment_name:
                writer.writerow(row)
                count += 1

print(f"Filtered {count} rows for experiment_name='{experiment_name}'", file=sys.stderr)
EOF
}

for ORIGIN in $ORIGINS; do
    ORIGIN_DIR="$WORK_DIR/$ORIGIN"
    mkdir -p "$ORIGIN_DIR"

    echo "============================================================"
    echo "Processing origin: $ORIGIN"
    echo "============================================================"

    # Paths for this origin
    ORIGIN_MGF="$MGF_DIR/${ORIGIN}.mgf"
    FILTERED_PREDS="$ORIGIN_DIR/${ORIGIN}_preds.csv"
    TARGET_FILE="$ORIGIN_DIR/${ORIGIN}_de_novo_results.csv"
    DB_FILE="$ORIGIN_DIR/${ORIGIN}_database_results.csv"

    if [ ! -f "$ORIGIN_MGF" ]; then
        echo "Warning: MGF file not found: $ORIGIN_MGF"
        echo "Skipping $ORIGIN"
        continue
    fi

    # Step 1: Filter combined predictions for the original (non-decoy) MGF
    echo "  Filtering predictions for $ORIGIN..."
    filter_csv "$COMBINED_PREDS" "$FILTERED_PREDS" "$ORIGIN"

    # Step 2: Preprocess target predictions
    echo "  Preprocessing target predictions..."
    novoboard preprocess \
        --denovo-file "$FILTERED_PREDS" \
        --denovo-output "$TARGET_FILE" \
        --db-mgf-file "$ORIGIN_MGF" \
        --db-output "$DB_FILE"

    # Step 3: Filter and preprocess decoy predictions
    DECOY_FILES=()
    LABELS=()
    for RATE in $SAMPLING_RATES; do
        DECOY_MGF="$MGF_DIR/${ORIGIN}_decoy_${RATE}.mgf"

        if [ ! -f "$DECOY_MGF" ]; then
            echo "  Warning: Decoy MGF not found: $DECOY_MGF"
            continue
        fi

        # The experiment_name for decoys should match the decoy MGF filename (without .mgf)
        DECOY_EXPERIMENT="${ORIGIN}_decoy_${RATE}"
        FILTERED_DECOY_PREDS="$ORIGIN_DIR/${DECOY_EXPERIMENT}_preds.csv"
        DECOY_DE_NOVO="$ORIGIN_DIR/${DECOY_EXPERIMENT}_de_novo_results.csv"

        echo "  Filtering predictions for ${DECOY_EXPERIMENT}..."
        filter_csv "$COMBINED_PREDS" "$FILTERED_DECOY_PREDS" "$DECOY_EXPERIMENT"

        # Check if we got any predictions
        DECOY_COUNT=$(wc -l < "$FILTERED_DECOY_PREDS")
        if [ "$DECOY_COUNT" -le 1 ]; then
            echo "  Warning: No predictions found for ${DECOY_EXPERIMENT}"
            continue
        fi

        echo "  Preprocessing decoy predictions..."
        novoboard preprocess \
            --denovo-file "$FILTERED_DECOY_PREDS" \
            --denovo-output "$DECOY_DE_NOVO" \
            --db-mgf-file "$DECOY_MGF" \
            --db-output /dev/null 2>/dev/null || true

        if [ -f "$DECOY_DE_NOVO" ]; then
            DECOY_FILES+=("$DECOY_DE_NOVO")
            # Convert rate to percentage label (e.g., 0.10 -> "10% Kept")
            # Note: RATE is the sampling/retention rate (peaks KEPT), not masking rate
            # decoy_0.10.mgf means 10% of peaks kept, 90% replaced with noise
            KEPT_PCT=$(echo "$RATE * 100" | bc | cut -d. -f1)
            LABELS+=("${KEPT_PCT}% Kept")
        fi
    done

    if [ ${#DECOY_FILES[@]} -eq 0 ]; then
        echo "Warning: No decoy files successfully processed for $ORIGIN"
        echo "Skipping $ORIGIN"
        continue
    fi

    OUTPUT_FILE="$OUTPUT_DIR/${ORIGIN}_fdr_validation.png"

    echo ""
    echo "  Target: $TARGET_FILE"
    echo "  Database: $DB_FILE"
    echo "  Decoys: ${#DECOY_FILES[@]} files"
    echo "  Output: $OUTPUT_FILE"
    echo ""

    # Step 4: Run novoboard fdr
    echo "  Running FDR analysis..."
    novoboard fdr \
        --target-file "$TARGET_FILE" \
        --decoy-files "${DECOY_FILES[@]}" \
        --labels "${LABELS[@]}" \
        --db-file "$DB_FILE" \
        --spectrum-file "$ORIGIN_MGF" \
        --output-file "$OUTPUT_FILE" \
        --score-column 'ALC (%)' \
        --aa-score-column 'local confidence (%)' \
        --fdr-max 0.1 \
        --dpi 300 \
        --tp-metric peptide

    echo ""
done

echo "============================================================"
echo "Done! Per-origin FDR plots saved to: $OUTPUT_DIR"
echo "Processed files saved to: $WORK_DIR"
echo "============================================================"
