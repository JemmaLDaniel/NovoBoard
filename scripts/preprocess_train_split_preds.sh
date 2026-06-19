#!/bin/bash
# Preprocess InstaNovo predictions for train-split datasets (celegans, sbrodae, PXD019483).
# For each origin and decoy rate, runs `novoboard preprocess` to produce:
#   - {prefix}_denovo.csv   (de novo results)
#   - {prefix}_db_output.csv (database ground truth from MGF labels)
#
# Usage:
#   bash scripts/preprocess_train_split_preds.sh 1000
#   bash scripts/preprocess_train_split_preds.sh 2000

set -e

SAMPLE_N="${1:?Usage: $0 <sample_n> (1000 or 2000)}"

TUNING_DIR="paper_analysis/tuning"
ORIGINS="celegans sbrodae PXD019483"
RATES="0.30 0.40 0.50 0.60 0.70"
PREFIX="annotated_train_sample${SAMPLE_N}"

for origin in $ORIGINS; do
    origin_dir="${TUNING_DIR}/${origin}"

    if [ ! -d "$origin_dir" ]; then
        echo "ERROR: directory not found: $origin_dir — skipping $origin"
        continue
    fi

    echo "============================================================"
    echo "  Origin: $origin"
    echo "============================================================"

    # --- Target (non-decoy) ---
    target_csv="${origin_dir}/${PREFIX}.csv"
    target_mgf="${origin_dir}/${PREFIX}.mgf"
    target_denovo="${origin_dir}/${PREFIX}_denovo.csv"
    target_db="${origin_dir}/${PREFIX}_db_output.csv"

    if [ ! -f "$target_csv" ]; then
        echo "  WARNING: target CSV not found: $target_csv — skipping target"
    elif [ ! -f "$target_mgf" ]; then
        echo "  WARNING: target MGF not found: $target_mgf — skipping target"
    else
        echo "  Preprocessing target predictions..."
        novoboard preprocess \
            --denovo-file "$target_csv" \
            --denovo-output "$target_denovo" \
            --db-mgf-file "$target_mgf" \
            --db-output "$target_db"
    fi

    # --- Decoys ---
    for rate in $RATES; do
        decoy_prefix="${PREFIX}_decoy_${rate}"
        decoy_csv="${origin_dir}/${decoy_prefix}.csv"
        decoy_mgf="${origin_dir}/${decoy_prefix}.mgf"
        decoy_denovo="${origin_dir}/${decoy_prefix}_denovo.csv"
        decoy_db="${origin_dir}/${decoy_prefix}_db_output.csv"

        if [ ! -f "$decoy_csv" ]; then
            echo "  WARNING: decoy CSV not found: $decoy_csv — skipping rate $rate"
            continue
        fi
        if [ ! -f "$decoy_mgf" ]; then
            echo "  WARNING: decoy MGF not found: $decoy_mgf — skipping rate $rate"
            continue
        fi

        echo "  Preprocessing decoy rate=${rate}..."
        novoboard preprocess \
            --denovo-file "$decoy_csv" \
            --denovo-output "$decoy_denovo" \
            --db-mgf-file "$decoy_mgf" \
            --db-output "$decoy_db"
    done

    echo ""
done

echo "Preprocessing complete!"
