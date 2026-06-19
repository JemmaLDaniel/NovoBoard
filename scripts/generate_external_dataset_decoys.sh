#!/usr/bin/env bash
# Generate MGF files and decoys for external datasets (celegans_raw, immuno2_raw).
#
# This script:
#   1. Exports each dataset from parquet to MGF format using sample_and_export_mgf.py
#      (celegans is sharded into N_SHARDS parts to avoid RAM OOMs during decoy generation)
#   2. Generates decoy spectra at RATE peaks kept using NovoBoard
#
# MGF title format: "{experiment_name} scan={scan_number}"  (matches per-origin pipeline)
# SCANS field:      0-based index within each shard file
#
# Usage:
#   bash scripts/generate_external_dataset_decoys.sh
#
# Environment variables:
#   CELEGANS_PARQUET  Path to celegans parquet (default: paper_datasets/celegans_raw.parquet)
#   IMMUNO2_PARQUET   Path to immuno2 parquet  (default: paper_datasets/immuno2_raw.parquet)
#   OUTPUT_DIR        Output directory for MGF files
#   RATE              Fraction of peaks to KEEP (default: 0.5  → 50% masked)
#   N_SHARDS          Number of shards for celegans MGF (default: 5)
#   SEED              Random seed (default: 42)

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CELEGANS_PARQUET="${CELEGANS_PARQUET:-/home/j-daniel/Documents/winnow/data/paper_datasets/celegans_raw.parquet}"
IMMUNO2_PARQUET="${IMMUNO2_PARQUET:-/home/j-daniel/Documents/winnow/data/paper_datasets/immuno2_raw.parquet}"
OUTPUT_DIR="${OUTPUT_DIR:-/home/j-daniel/Documents/winnow/novoboard/benchmark_data/instanovo_inputs/external_datasets}"
RATE="${RATE:-0.5}"
N_SHARDS="${N_SHARDS:-5}"
SEED="${SEED:-42}"

mkdir -p "${OUTPUT_DIR}"

echo "=============================================="
echo "External Dataset Decoy Generation Pipeline"
echo "=============================================="
echo "celegans parquet : ${CELEGANS_PARQUET}"
echo "immuno2  parquet : ${IMMUNO2_PARQUET}"
echo "Output directory : ${OUTPUT_DIR}"
echo "Sampling rate    : ${RATE} peaks kept ($(echo "scale=0; (1-${RATE})*100/1" | bc)% masked)"
echo "celegans shards  : ${N_SHARDS}"
echo "Random seed      : ${SEED}"
echo "=============================================="
echo ""

# ---------------------------------------------------------------------------
# celegans_raw  (sharded to keep each novoboard decoy run within RAM budget)
# ---------------------------------------------------------------------------
echo "======================================"
echo "Processing: celegans_raw (${N_SHARDS} shards)"
echo "======================================"

for shard_idx in $(seq 0 $((N_SHARDS - 1))); do
    shard_num=$((shard_idx + 1))
    MGF_FILE="${OUTPUT_DIR}/celegans_raw_${shard_num}_${N_SHARDS}.mgf"

    echo ""
    echo "  [celegans shard ${shard_num}/${N_SHARDS}] Exporting MGF..."
    uv run python scripts/sample_and_export_mgf.py \
        --parquet "${CELEGANS_PARQUET}" \
        --output  "${MGF_FILE}" \
        --n-shards "${N_SHARDS}" \
        --shard-idx "${shard_idx}" \
        --seed "${SEED}"

    echo "  [celegans shard ${shard_num}/${N_SHARDS}] Generating decoys (rate=${RATE})..."
    uv run novoboard decoy \
        --spectrum-file "${MGF_FILE}" \
        --sampling-strategy random \
        --sampling-rate "${RATE}" \
        --seed "${SEED}"

    echo "  [celegans shard ${shard_num}/${N_SHARDS}] Done."
done

echo ""
echo "======================================"
echo "Processing: immuno2_raw"
echo "======================================"

IMMUNO2_MGF="${OUTPUT_DIR}/immuno2_raw.mgf"

echo ""
echo "  [immuno2] Exporting MGF..."
uv run python scripts/sample_and_export_mgf.py \
    --parquet "${IMMUNO2_PARQUET}" \
    --output  "${IMMUNO2_MGF}" \
    --seed "${SEED}"

echo "  [immuno2] Generating decoys (rate=${RATE})..."
uv run novoboard decoy \
    --spectrum-file "${IMMUNO2_MGF}" \
    --sampling-strategy random \
    --sampling-rate "${RATE}" \
    --seed "${SEED}"

echo "  [immuno2] Done."

echo ""
echo "=============================================="
echo "Pipeline complete!"
echo "=============================================="
echo "Output files in: ${OUTPUT_DIR}"
echo ""
echo "Generated files:"
echo "  celegans_raw_{1..${N_SHARDS}}_${N_SHARDS}.mgf          (target spectra)"
echo "  celegans_raw_{1..${N_SHARDS}}_${N_SHARDS}_decoy_${RATE}.mgf (decoys)"
echo "  immuno2_raw.mgf                      (target spectra)"
echo "  immuno2_raw_decoy_${RATE}.mgf         (decoys)"
