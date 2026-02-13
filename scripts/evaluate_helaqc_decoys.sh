#!/bin/bash

novoboard fdr \
    --target-file data/converted/helaqc_de_novo_results.csv \
    --decoy-files \
        data/converted/helaqc_decoy_0.10_de_novo_results.csv \
        data/converted/helaqc_decoy_0.20_de_novo_results.csv \
        data/converted/helaqc_decoy_0.30_de_novo_results.csv \
        data/converted/helaqc_decoy_0.40_de_novo_results.csv \
        data/converted/helaqc_decoy_0.50_de_novo_results.csv \
        data/converted/helaqc_decoy_0.60_de_novo_results.csv \
        data/converted/helaqc_decoy_0.70_de_novo_results.csv \
        data/converted/helaqc_decoy_0.80_de_novo_results.csv \
        data/converted/helaqc_decoy_0.90_de_novo_results.csv \
    --labels \
        "10% Kept" "20% Kept" "30% Kept" "40% Kept" "50% Kept" \
        "60% Kept" "70% Kept" "80% Kept" "90% Kept" \
    --db-file data/converted/helaqc_database_results.csv \
    --spectrum-file /home/j-daniel/Documents/winnow/novoboard/helaqc.mgf \
    --output-file fig/helaqc_fdr_validation.png \
    --score-column "ALC (%)" \
    --aa-score-column "local confidence (%)" \
    --fdr-max 0.1 \
    --dpi 300 \
    --tp-metric peptide \
    --no-monotonic
