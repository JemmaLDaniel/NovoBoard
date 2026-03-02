#!/usr/bin/env bash
# Generate per-origin MGF files and decoys at multiple masking rates.
#
# This script:
# 1. Samples spectra from each origin in a parquet file
# 2. Exports them to MGF format
# 3. Generates decoy spectra at multiple masking rates using NovoBoard
#
# Usage:
#   bash scripts/generate_per_origin_decoys.sh
#
# Environment variables:
#   PARQUET_PATH: Path to input parquet file (default: train.parquet)
#   OUTPUT_DIR: Directory for output MGF files (default: benchmark_data/instanovo_inputs)
#   SAMPLE_N: Number of spectra to sample per origin (default: 2000)
#   SEED: Random seed for sampling (default: 42)

set -euo pipefail

# Configuration
# PARQUET_PATH can be a single file, multiple files (space-separated), or a glob pattern
PARQUET_PATH="${PARQUET_PATH:-/home/j-daniel/Documents/winnow/data/paper_datasets/train.parquet /home/j-daniel/Documents/winnow/data/paper_datasets/val.parquet /home/j-daniel/Documents/winnow/data/paper_datasets/test.parquet}"
OUTPUT_DIR="${OUTPUT_DIR:-/home/j-daniel/Documents/winnow/novoboard/benchmark_data/per_origin_mgf}"
SAMPLE_N="${SAMPLE_N:-2000}"
SEED="${SEED:-42}"

# Origins to process (all origins in train.parquet)
ORIGINS=(
    "hepg2"
    "gluc"
    "tplantibodies"
    "helaqc"
    "sbrodae"
    "snakevenoms"
    "woundfluids"
    "immuno"
    "herceptin"
)

# Masking rates to generate (fraction of peaks to KEEP)
# 0.1 = 90% masked (10% kept), 0.9 = 10% masked (90% kept)
# Can be overridden via RATES environment variable
if [[ -n "${RATES:-}" ]]; then
    # Convert space-separated string to array
    read -ra RATES_ARRAY <<< "${RATES}"
else
    RATES_ARRAY=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Count parquet files
# shellcheck disable=SC2086
PARQUET_COUNT=$(echo ${PARQUET_PATH} | wc -w)

echo "=============================================="
echo "Per-Origin Decoy Generation Pipeline"
echo "=============================================="
echo "Input parquet files (${PARQUET_COUNT}):"
for pq in ${PARQUET_PATH}; do
    echo "  - ${pq}"
done
echo "Output directory: ${OUTPUT_DIR}"
echo "Sample size per origin: ${SAMPLE_N}"
echo "Random seed: ${SEED}"
echo "Origins: ${ORIGINS[*]}"
echo "Masking rates (peaks kept): ${RATES_ARRAY[*]}"
echo "=============================================="
echo ""

# Process each origin
for origin in "${ORIGINS[@]}"; do
    echo ""
    echo "======================================"
    echo "Processing origin: ${origin}"
    echo "======================================"

    MGF_FILE="${OUTPUT_DIR}/${origin}.mgf"

    # Step 1: Export MGF from parquet (PARQUET_PATH may contain multiple files)
    echo "Step 1: Exporting MGF..."
    # shellcheck disable=SC2086
    python scripts/sample_and_export_mgf.py \
        --parquet ${PARQUET_PATH} \
        --output "${MGF_FILE}" \
        --origin "${origin}" \
        --sample-n "${SAMPLE_N}" \
        --seed "${SEED}"

    # Step 2: Generate decoys at each masking rate
    echo "Step 2: Generating decoys..."
    for rate in "${RATES_ARRAY[@]}"; do
        echo "  Generating decoys at ${rate} peaks kept..."
        novoboard decoy \
            --spectrum-file "${MGF_FILE}" \
            --sampling-strategy random \
            --sampling-rate "${rate}" \
            --seed "${SEED}"
    done

    echo "Completed ${origin}"
done

echo ""
echo "=============================================="
echo "Pipeline complete!"
echo "=============================================="
echo "Output files in: ${OUTPUT_DIR}"
echo ""
echo "Generated files per origin:"
echo "  - {origin}.mgf (sampled spectra)"
echo "  - {origin}_decoy_{rate}.mgf (decoys at each rate)"
