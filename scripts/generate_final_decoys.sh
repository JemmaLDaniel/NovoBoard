#!/bin/bash
# Generate decoys for the full winnow train/val/unlabelled MGF splits.
#
# For each origin, generates decoy MGFs at a per-origin masking rate for:
#   - annotated_train.mgf
#   - annotated_val.mgf
#   - raw_unlabelled.mgf
#
# Masking rate = fraction of peaks REMOVED (replaced with noise).
# Internally converted to kept rate for novoboard decoy --sampling-rate.
#
# Output goes to {WINNOW_DIR}/{origin}_split_mgf_final_decoys/
#
# Usage:
#   bash scripts/generate_final_decoys.sh celegans=0.40 sbrodae=0.50 PXD019483=0.40

set -e

if [ $# -eq 0 ]; then
    echo "Usage: $0 origin1=masking_rate1 origin2=masking_rate2 ..."
    echo "  masking_rate = fraction of peaks REMOVED (e.g. 0.40 = 40% masked)"
    echo "Example: $0 celegans=0.40 sbrodae=0.50 PXD019483=0.40"
    exit 1
fi

WINNOW_DIR="/home/j-daniel/repos/winnow"
MGFS="annotated_test.mgf annotated_val.mgf raw_unlabelled.mgf"
SEED=42

for arg in "$@"; do
    origin="${arg%%=*}"
    masking_rate="${arg##*=}"

    if [ "$origin" = "$masking_rate" ]; then
        echo "ERROR: invalid argument '$arg' — expected origin=masking_rate (e.g. celegans=0.40)"
        exit 1
    fi

    kept_rate=$(echo "1.0 - $masking_rate" | bc)

    src_dir="${WINNOW_DIR}/${origin}_split_mgf"
    dst_dir="${WINNOW_DIR}/${origin}_split_mgf_final_decoys"

    echo "============================================================"
    echo "  Origin: $origin  (masking=${masking_rate}, kept=${kept_rate})"
    echo "============================================================"

    mkdir -p "$dst_dir"

    for mgf_name in $MGFS; do
        src_mgf="${src_dir}/${mgf_name}"

        if [ ! -f "$src_mgf" ]; then
            echo "  WARNING: not found: $src_mgf — skipping"
            continue
        fi

        dst_mgf="${dst_dir}/${mgf_name}"
        if [ ! -f "$dst_mgf" ]; then
            echo "  Copying $mgf_name to output dir..."
            cp "$src_mgf" "$dst_mgf"
        fi

        echo "  Generating decoy for $mgf_name (masking=${masking_rate})..."
        novoboard decoy \
            --spectrum-file "$dst_mgf" \
            --sampling-strategy random \
            --sampling-rate "$kept_rate" \
            --seed $SEED
    done

    echo ""
done

echo "Done!"
