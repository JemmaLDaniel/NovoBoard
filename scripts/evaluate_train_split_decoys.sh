#!/bin/bash
# Evaluate FDR for train-split datasets (celegans, sbrodae, PXD019483).
#
# For each origin, runs `novoboard fdr` comparing decoy-based FDR estimates
# against database-grounded truth at multiple masking rates, producing a
# validation plot per origin.
#
# Prerequisites:
#   Run `scripts/preprocess_train_split_preds.sh <sample_n>` first to generate
#   *_denovo.csv and *_db_output.csv files.
#
# Usage:
#   bash scripts/evaluate_train_split_decoys.sh 1000
#   bash scripts/evaluate_train_split_decoys.sh 2000

set -e

SAMPLE_N="${1:?Usage: $0 <sample_n> (1000 or 2000)}"

TUNING_DIR="paper_analysis/tuning"
ORIGINS="celegans sbrodae PXD019483"
RATES="0.30 0.40 0.50 0.60 0.70"
PREFIX="annotated_train_sample${SAMPLE_N}"
OUTPUT_DIR="paper_analysis/fig/tuning_n${SAMPLE_N}"

mkdir -p "$OUTPUT_DIR"

for origin in $ORIGINS; do
    origin_dir="${TUNING_DIR}/${origin}"

    echo "============================================================"
    echo "  Origin: $origin"
    echo "============================================================"

    target_denovo="${origin_dir}/${PREFIX}_denovo.csv"
    target_db="${origin_dir}/${PREFIX}_db_output.csv"
    spectrum_file="${origin_dir}/${PREFIX}.mgf"

    if [ ! -f "$target_denovo" ]; then
        echo "  ERROR: target denovo not found: $target_denovo — run preprocess first"
        continue
    fi
    if [ ! -f "$target_db" ]; then
        echo "  ERROR: target db not found: $target_db — run preprocess first"
        continue
    fi
    if [ ! -f "$spectrum_file" ]; then
        echo "  ERROR: spectrum file not found: $spectrum_file"
        continue
    fi

    DECOY_FILES=()
    LABELS=()
    for rate in $RATES; do
        decoy_denovo="${origin_dir}/${PREFIX}_decoy_${rate}_denovo.csv"
        if [ ! -f "$decoy_denovo" ]; then
            echo "  WARNING: decoy denovo not found: $decoy_denovo — skipping rate $rate"
            continue
        fi
        DECOY_FILES+=("$decoy_denovo")
        kept_pct=$(echo "$rate * 100" | bc | cut -d. -f1)
        masked_pct=$((100 - kept_pct))
        LABELS+=("${masked_pct}% Masked")
    done

    if [ ${#DECOY_FILES[@]} -eq 0 ]; then
        echo "  WARNING: no decoy files found for $origin — skipping"
        continue
    fi

    output_file="${OUTPUT_DIR}/${origin}_fdr_validation.png"

    echo "  Target:   $target_denovo"
    echo "  Database: $target_db"
    echo "  Spectrum: $spectrum_file"
    echo "  Decoys:   ${#DECOY_FILES[@]} files"
    echo "  Output:   $output_file"
    echo ""

    novoboard fdr \
        --target-file "$target_denovo" \
        --decoy-files "${DECOY_FILES[@]}" \
        --labels "${LABELS[@]}" \
        --db-file "$target_db" \
        --spectrum-file "$spectrum_file" \
        --output-file "$output_file" \
        --score-column "ALC (%)" \
        --aa-score-column "local confidence (%)" \
        --fdr-max 0.1 \
        --dpi 300 \
        --tp-metric peptide

    echo ""
done

echo "============================================================"
echo "Done! FDR validation plots saved to: $OUTPUT_DIR"
echo "============================================================"
