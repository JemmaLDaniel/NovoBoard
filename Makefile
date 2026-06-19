# NovoBoard Makefile
# Commands for running common workflows

.PHONY: help test fdr-combined fdr-per-origin separate-by-origin evaluate-per-origin \
        extract-metrics visualize-masking clean-per-origin clean-cache \
        per-origin-decoys single-origin-decoys export-mgf train-split-decoys \
        fdr-generalizability-analysis spectrum-distributions confidence-distributions \
        peptide-characteristics external-dataset-decoys \
        preprocess-train-splits fdr-train-splits select-masking-rate \
        final-decoys \
        preprocess-datasets fdr-datasets \
        $(foreach o,$(DATASET_ORIGINS),preprocess-dataset-$(o) fdr-dataset-$(o))

# Default target
help:
	@echo "NovoBoard Makefile"
	@echo ""
	@echo "Available targets:"
	@echo ""
	@echo "  Testing:"
	@echo "    test              - Run pytest test suite"
	@echo ""
	@echo "  FDR Analysis:"
	@echo "    fdr-generalizability-analysis - Run the FULL generalizability analysis pipeline"
	@echo "                          (decoys → FDR eval → metrics → figures → characteristics)"
	@echo "    fdr-combined      - Run FDR analysis on combined train_sample10 data"
	@echo "    fdr-per-origin    - Run complete per-origin FDR analysis pipeline"
	@echo "    separate-by-origin - Separate train_sample10 data by origin"
	@echo "    evaluate-per-origin - Evaluate FDR for each origin"
	@echo "    extract-metrics   - Extract FDR metrics (MAE/Bias/Slope/R²/True@X%)"
	@echo "    spectrum-distributions - Generate intensity/m/z distribution figures"
	@echo "                          Usage: make spectrum-distributions VERIFY_ORIGINS='hepg2 gluc helaqc'"
	@echo "    confidence-distributions - Extract TP/FP/Decoy confidence score tables"
	@echo "                          Usage: make confidence-distributions CONFIDENCE_RATES='0.50 0.81'"
	@echo "    peptide-characteristics - Extract peptide length & charge distribution table"
	@echo "    external-dataset-decoys - Generate MGF files and decoys for external datasets"
	@echo "                          (celegans sharded; immuno2 whole)"
	@echo "                          Usage: make external-dataset-decoys RATE=0.4 N_SHARDS=5"
	@echo "    train-split-decoys  - Sample & generate decoys for winnow split MGFs"
	@echo "                          (sbrodae, celegans, PXD019483)"
	@echo "                          Usage: make train-split-decoys TRAIN_SPLIT_N=1000 TRAIN_SPLIT_RATES='0.3 0.5 0.7'"
	@echo "    preprocess-train-splits - Preprocess InstaNovo preds for train-split tuning"
	@echo "                          Usage: make preprocess-train-splits TRAIN_SPLIT_SAMPLE_N=2000"
	@echo "    fdr-train-splits  - Preprocess + evaluate FDR for train-split tuning"
	@echo "                          Usage: make fdr-train-splits TRAIN_SPLIT_SAMPLE_N=2000"
	@echo "    select-masking-rate - Find best masking rate per origin via MSE"
	@echo "                          Usage: make select-masking-rate FDR_MAX=0.1"
	@echo "    final-decoys      - Generate decoys for full winnow train/val/unlabelled splits"
	@echo "                          Usage: make final-decoys FINAL_MASKING_RATES='celegans=0.30 sbrodae=0.50 PXD019483=0.30'"
	@echo ""
	@echo "  Dataset FDR Estimation (de novo, no database):"
	@echo "    preprocess-datasets - Preprocess InstaNovo CSVs in datasets/ to NovoBoard format"
	@echo "    fdr-datasets        - Estimate per-PSM FDR & q-values for annotated_test + raw_unlabelled"
	@echo ""
	@echo "  Visualization:"
	@echo "    visualize-masking - Compare original spectra with pre-generated decoys"
	@echo "                        Usage: make visualize-masking ORIGIN=gluc DECOY_RATES='0.10 0.50'"
	@echo ""
	@echo "  Decoy Generation:"
	@echo "    per-origin-decoys   - Generate per-origin MGF and decoy files"
	@echo "                          Usage: make per-origin-decoys RATES='0.3 0.5 0.7'"
	@echo "    single-origin-decoys - Generate MGF and decoys for one origin at one rate"
	@echo "                          Usage: make single-origin-decoys ORIGIN=hepg2 RATE=0.7"
	@echo "    export-mgf          - Export MGF from parquet (no decoys)"
	@echo ""
	@echo "  Cleanup:"
	@echo "    clean-cache       - Clear accuracy cache files (forces recalculation)"
	@echo "    clean-per-origin  - Remove per-origin output files"
	@echo ""
	@echo "Example usage:"
	@echo "  make test"
	@echo "  make fdr-combined"
	@echo "  make fdr-per-origin"
	@echo "  make visualize-masking ORIGIN=gluc DECOY_RATES='0.10 0.50' N=3"
	@echo "  make visualize-masking ORIGIN=hepg2 DECOY_RATES='0.30' N=5"

# =============================================================================
# Testing
# =============================================================================

test:
	@echo "Running pytest..."
	uv run pytest tests/ -v

# =============================================================================
# FDR Analysis
# =============================================================================

# Run FDR analysis on combined train_sample10 data
fdr-combined:
	@echo "Running FDR analysis on combined train_sample10 data..."
	uv run bash scripts/evaluate_train_decoys.sh

# Complete per-origin FDR analysis pipeline
fdr-per-origin: separate-by-origin evaluate-per-origin
	@echo "Per-origin FDR analysis complete!"

# Step 1: Separate data by origin
separate-by-origin:
	@echo "Separating train_sample10 data by origin..."
	uv run python scripts/separate_by_origin.py

# Step 2: Evaluate FDR for each origin
evaluate-per-origin:
	@echo "Evaluating FDR for each origin..."
	uv run bash scripts/evaluate_train_decoys_per_origin.sh

# Extract FDR metrics for generalizability analysis
extract-metrics:
	@echo "Extracting FDR metrics..."
	uv run python scripts/extract_fdr_metrics.py --print-tables
	@echo "Metrics saved to analysis/"

# Verify decoy spectrum quality and generate intensity/m/z distribution figures
# Covers: fig/spectrum_distributions/ (Section 6.2.2 of fdr_generalizability_analysis.md)
VERIFY_ORIGINS ?= hepg2 gluc helaqc
spectrum-distributions:
	@echo "Generating spectrum distribution figures for: $(VERIFY_ORIGINS)..."
	@for origin in $(VERIFY_ORIGINS); do \
		echo "  Processing $$origin..."; \
		uv run python scripts/verify_decoy_spectra.py --origin $$origin --fig-output-dir fig/spectrum_distributions --csv-output-dir analysis; \
	done
	@echo "Spectrum distribution figures saved to fig/spectrum_distributions/"

# Extract confidence score distributions (TP / FP / decoy) per origin
# Covers: Table in Section 6.2.4 of fdr_generalizability_analysis.md
CONFIDENCE_RATES ?= 0.50 0.80
confidence-distributions:
	@echo "Extracting confidence score distributions..."
	@for rate in $(CONFIDENCE_RATES); do \
		echo "  Rate $$rate (kept)..."; \
		uv run python scripts/extract_confidence_distributions.py --retention-rate $$rate; \
	done
	@echo "Confidence distribution tables saved to analysis/"

# Extract peptide characteristics per origin (length, charge distribution)
# Covers: Table in Section 6.2.3 of fdr_generalizability_analysis.md
PAPER_PARQUET ?= /home/j-daniel/Documents/winnow/data/paper_datasets/train.parquet \
                 /home/j-daniel/Documents/winnow/data/paper_datasets/val.parquet \
                 /home/j-daniel/Documents/winnow/data/paper_datasets/test.parquet
peptide-characteristics:
	@echo "Extracting peptide characteristics per origin..."
	uv run python scripts/extract_peptide_characteristics.py \
		--parquet-files $(PAPER_PARQUET) \
		--print-tables
	@echo "Peptide characteristics saved to analysis/"

# Generate MGF files and decoys for external datasets (celegans_raw, immuno2_raw)
# celegans is sharded into N_SHARDS parts to keep novoboard decoy RAM usage bounded.
# Usage:
#   make external-dataset-decoys                     # defaults: RATE=0.4, N_SHARDS=5
#   make external-dataset-decoys RATE=0.5 N_SHARDS=10
EXTERNAL_OUTPUT_DIR ?= /home/j-daniel/Documents/winnow/novoboard/benchmark_data/instanovo_inputs/external_datasets
CELEGANS_PARQUET ?= /home/j-daniel/Documents/winnow/data/paper_datasets/celegans_raw.parquet
IMMUNO2_PARQUET  ?= /home/j-daniel/Documents/winnow/data/paper_datasets/immuno2_raw.parquet
EXTERNAL_RATE    ?= 0.5
N_SHARDS         ?= 5

external-dataset-decoys:
	@echo "Generating MGF files and decoys for external datasets (rate=$(EXTERNAL_RATE))..."
	CELEGANS_PARQUET=$(CELEGANS_PARQUET) \
	IMMUNO2_PARQUET=$(IMMUNO2_PARQUET) \
	OUTPUT_DIR=$(EXTERNAL_OUTPUT_DIR) \
	RATE=$(EXTERNAL_RATE) \
	N_SHARDS=$(N_SHARDS) \
	SEED=$(SEED) \
		uv run bash scripts/generate_external_dataset_decoys.sh

# Run the full generalizability analysis pipeline
# Reproduces all tables and figures in analysis/fdr_generalizability_analysis.md
#
# Steps:
#   1. fdr-per-origin            - Separate combined predictions by origin,
#    								then run FDR evaluation for each origin/rate.
#   2. extract-metrics           - Compute MAE/Bias/Slope/R²/True@X% tables (Tables 2-6)
#   3. spectrum-distributions    - Intensity/m/z distribution figures (Section 6.2.2)
#   4. confidence-distributions  - TP/FP/Decoy confidence tables (Section 6.2.4)
#   5. peptide-characteristics   - Peptide length & charge distribution table (Section 6.2.3)
#   6. visualize-masking         - Example masking figures for gluc and hepg2 (Section 1 & 6.1)
#
# Usage:
#   make fdr-generalizability-analysis
#   make fdr-generalizability-analysis RATES="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9" SAMPLE_N=2000
fdr-generalizability-analysis: fdr-per-origin extract-metrics \
                                spectrum-distributions confidence-distributions \
                                peptide-characteristics
	@echo "Generating example masking visualizations..."
	$(MAKE) visualize-masking ORIGIN=gluc DECOY_RATES="0.30 0.50 0.70" N=3 SEED=42
	$(MAKE) visualize-masking ORIGIN=hepg2 DECOY_RATES="0.30 0.50 0.70" N=3 SEED=42
	@echo ""
	@echo "=============================================="
	@echo "Generalizability analysis complete!"
	@echo "Outputs:"
	@echo "  Tables  -> analysis/"
	@echo "  Figures -> fig/per_origin/, fig/spectrum_distributions/, fig/masking_examples/"
	@echo "=============================================="

# Preprocess InstaNovo predictions for train-split tuning datasets
# (celegans, sbrodae, PXD019483) at all decoy rates
# Usage: make preprocess-train-splits TRAIN_SPLIT_SAMPLE_N=2000
TRAIN_SPLIT_SAMPLE_N ?= 1000
preprocess-train-splits:
	@echo "Preprocessing train-split predictions (n=$(TRAIN_SPLIT_SAMPLE_N))..."
	uv run bash scripts/preprocess_train_split_preds.sh $(TRAIN_SPLIT_SAMPLE_N)

# Preprocess + evaluate FDR for train-split tuning datasets
# Usage: make fdr-train-splits TRAIN_SPLIT_SAMPLE_N=2000
fdr-train-splits: preprocess-train-splits
	@echo "Evaluating FDR for train-split datasets (n=$(TRAIN_SPLIT_SAMPLE_N))..."
	uv run bash scripts/evaluate_train_split_decoys.sh $(TRAIN_SPLIT_SAMPLE_N)

# Select best masking rate per origin based on MSE from ideal calibration
# Usage: make select-masking-rate TRAIN_SPLIT_SAMPLE_N=1000
#        make select-masking-rate TRAIN_SPLIT_SAMPLE_N=2000 FDR_MAX=0.1
FDR_MAX ?=
select-masking-rate:
	uv run python scripts/select_masking_rate.py $(TRAIN_SPLIT_SAMPLE_N) \
		$(if $(FDR_MAX),--fdr-max $(FDR_MAX))

# Generate final decoys for full winnow train/val/unlabelled splits
# Copies source MGFs into {origin}_split_mgf_final_decoys/ and generates
# decoy MGFs alongside them at per-origin masking rates (fraction REMOVED).
# Usage: make final-decoys FINAL_MASKING_RATES="celegans=0.30 sbrodae=0.50 PXD019483=0.30"
FINAL_MASKING_RATES ?= celegans=0.30 sbrodae=0.50 PXD019483=0.30 helaqc=0.50
final-decoys:
	@echo "Generating final decoys..."
	uv run bash scripts/generate_final_decoys.sh $(FINAL_MASKING_RATES)

# =============================================================================
# Visualization
# =============================================================================

# Default values for visualization
ORIGIN ?= hepg2
N ?= 3
SEED ?= 42
MGF_DIR ?= /home/j-daniel/Documents/winnow/novoboard/benchmark_data/instanovo_inputs/per_origin_mgf
# DECOY_RATES: space-separated list of rates to visualize (fraction of peaks KEPT)
DECOY_RATES ?= 0.10 0.50

# Visualize original spectra alongside pre-generated decoy spectra
# Usage: make visualize-masking ORIGIN=gluc DECOY_RATES="0.10 0.50" N=3
visualize-masking:
	@echo "Visualizing $(N) spectra from $(ORIGIN) with decoys at rates: $(DECOY_RATES)"
	$(eval MGF_FILE := $(MGF_DIR)/$(ORIGIN).mgf)
	$(eval DECOY_FILES := $(foreach rate,$(DECOY_RATES),$(MGF_DIR)/$(ORIGIN)_decoy_$(rate).mgf))
	uv run python scripts/visualize_masking.py \
		--mgf-file $(MGF_FILE) \
		--decoy-files $(DECOY_FILES) \
		--n-spectra $(N) \
		--seed $(SEED)
	@echo "Plots saved to fig/masking_examples/"

# =============================================================================
# Cleanup
# =============================================================================

# Clean up per-origin output files
clean-per-origin:
	@echo "Removing per-origin output files..."
	rm -rf /home/j-daniel/Documents/winnow/novoboard/benchmark_data/per_origin
	rm -rf fig/per_origin
	@echo "Cleaned per-origin outputs."

# Clear accuracy cache files (forces recalculation on next FDR run)
clean-cache:
	@echo "Clearing accuracy cache files..."
	rm -f /home/j-daniel/Documents/winnow/novoboard/benchmark_data/novoboard_inputs/*_accuracy.csv
	rm -f /home/j-daniel/Documents/winnow/novoboard/benchmark_data/per_origin/*/*_accuracy.csv
	@echo "Accuracy cache cleared. Next FDR run will recalculate from MGF files."

# =============================================================================
# Decoy Generation
# =============================================================================

# Default values for decoy generation (space-separated list of parquet files)
PARQUET ?= /home/j-daniel/Documents/winnow/data/paper_datasets/train.parquet /home/j-daniel/Documents/winnow/data/paper_datasets/val.parquet /home/j-daniel/Documents/winnow/data/paper_datasets/test.parquet
OUTPUT_DIR ?= /home/j-daniel/Documents/winnow/novoboard/benchmark_data/per_origin_mgf
SAMPLE_N ?= 2000
# Retention rates (fraction of peaks to KEEP): space-separated list
RATES ?= 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
# Single retention rate for single-origin-decoys
RATE ?= 0.5

# Generate per-origin MGF files and decoys for all origins
# This samples from train/val/test.parquet and generates decoys at specified masking rates
# Usage: make per-origin-decoys RATES="0.3 0.5 0.7"
per-origin-decoys:
	@echo "Generating per-origin MGF files and decoys..."
	PARQUET_PATH="$(PARQUET)" OUTPUT_DIR=$(OUTPUT_DIR) SAMPLE_N=$(SAMPLE_N) SEED=$(SEED) RATES="$(RATES)" \
		uv run bash scripts/generate_per_origin_decoys.sh

# Generate MGF and decoys for a single origin at a single retention rate
# Usage: make single-origin-decoys ORIGIN=hepg2 RATE=0.7 SAMPLE_N=2000
single-origin-decoys:
	@echo "Generating decoys for $(ORIGIN) at rate $(RATE)..."
	@mkdir -p $(OUTPUT_DIR)
	uv run python scripts/sample_and_export_mgf.py \
		--parquet $(PARQUET) \
		--output $(OUTPUT_DIR)/$(ORIGIN).mgf \
		--origin $(ORIGIN) \
		--sample-n $(SAMPLE_N) \
		--seed $(SEED)
	uv run novoboard decoy \
		--spectrum-file $(OUTPUT_DIR)/$(ORIGIN).mgf \
		--sampling-strategy random \
		--sampling-rate $(RATE) \
		--seed $(SEED)

# Export MGF from parquet file (without generating decoys)
# Usage: make export-mgf PARQUET=/path/to/data.parquet OUTPUT=/path/to/output.mgf ORIGIN=hepg2 SAMPLE_N=2000
OUTPUT ?= output.mgf
export-mgf:
	@echo "Exporting MGF from $(PARQUET) to $(OUTPUT)..."
	uv run python scripts/sample_and_export_mgf.py \
		--parquet $(PARQUET) \
		--output $(OUTPUT) \
		$(if $(ORIGIN),--origin $(ORIGIN)) \
		$(if $(SAMPLE_N),--sample-n $(SAMPLE_N)) \
		--seed $(SEED)

# Generate decoys for pre-split annotated_train.mgf files from winnow
# Samples TRAIN_SPLIT_N spectra from each annotated_train.mgf, saves as
# annotated_train_sample{N}.mgf in {origin}_split_mgf_decoys/, then generates
# matching decoy files alongside it.
# Usage:
#   make train-split-decoys                        # defaults: TRAIN_SPLIT_RATE=0.5, N=1000
#   make train-split-decoys TRAIN_SPLIT_RATE=0.4
#   make train-split-decoys TRAIN_SPLIT_RATES="0.3 0.5 0.7"
#   make train-split-decoys TRAIN_SPLIT_N=500
WINNOW_SPLIT_DIR    ?= /home/j-daniel/repos/winnow
TRAIN_SPLIT_ORIGINS ?= sbrodae celegans PXD019483
TRAIN_SPLIT_N       ?= 1000
TRAIN_SPLIT_RATE    ?= 0.5
TRAIN_SPLIT_RATES   ?= $(TRAIN_SPLIT_RATE)

train-split-decoys:
	@for origin in $(TRAIN_SPLIT_ORIGINS); do \
		src_dir="$(WINNOW_SPLIT_DIR)/$${origin}_split_mgf"; \
		dst_dir="$(WINNOW_SPLIT_DIR)/$${origin}_split_mgf_decoys"; \
		src_mgf="$${src_dir}/annotated_train.mgf"; \
		sampled_mgf="$${dst_dir}/annotated_train_sample$(TRAIN_SPLIT_N).mgf"; \
		if [ ! -f "$$src_mgf" ]; then \
			echo "ERROR: $$src_mgf not found, skipping $$origin"; \
			continue; \
		fi; \
		echo "=== $$origin ==="; \
		mkdir -p "$$dst_dir"; \
		echo "  Sampling $(TRAIN_SPLIT_N) spectra..."; \
		uv run python scripts/sample_mgf.py \
			--input "$$src_mgf" \
			--output "$$sampled_mgf" \
			--n $(TRAIN_SPLIT_N) \
			--seed $(SEED); \
		for rate in $(TRAIN_SPLIT_RATES); do \
			echo "  Generating decoy at rate $$rate..."; \
			uv run novoboard decoy \
				--spectrum-file "$$sampled_mgf" \
				--sampling-strategy random \
				--sampling-rate $$rate \
				--seed $(SEED); \
		done; \
	done
	@echo "Train split decoys complete!"

# =============================================================================
# Dataset FDR Estimation (de novo, no database validation)
# =============================================================================

# Per-origin decoy rates (fraction of peaks KEPT in the decoy spectrum filename)
DATASET_DIR         ?= datasets
DATASET_ORIGINS     ?= celegans sbrodae PXD019483 helaqc
DATASET_SPLITS      ?= annotated_test raw_unlabelled
DATASET_NOVOBOARD   ?= novoboard

# Map origin -> kept rate used in decoy filenames
DECOY_RATE_celegans   ?= 0.70
DECOY_RATE_sbrodae    ?= 0.50
DECOY_RATE_PXD019483  ?= 0.70
DECOY_RATE_helaqc     ?= 0.50

# Preprocess all raw InstaNovo CSVs in datasets/ to NovoBoard format.
# For each origin and split, converts {split}.csv and {split}_decoy_{rate}.csv
# into datasets/{origin}/novoboard/{split}.csv etc.
#
# Uses Make-level $(eval) to resolve the per-origin DECOY_RATE_<origin> variable,
# avoiding bash-only ${!var} indirect expansion.
define preprocess_origin
# $(1) = origin name
preprocess-datasets:: preprocess-dataset-$(1)
preprocess-dataset-$(1):
	@echo "=== Preprocessing $(1) ==="
	@mkdir -p "$(DATASET_DIR)/$(1)/$(DATASET_NOVOBOARD)"
	@for split in $(DATASET_SPLITS); do \
		src="$(DATASET_DIR)/$(1)/$$$${split}.csv"; \
		decoy_src="$(DATASET_DIR)/$(1)/$$$${split}_decoy_$(DECOY_RATE_$(1)).csv"; \
		dst="$(DATASET_DIR)/$(1)/$(DATASET_NOVOBOARD)/$$$${split}.csv"; \
		decoy_dst="$(DATASET_DIR)/$(1)/$(DATASET_NOVOBOARD)/$$$${split}_decoy_$(DECOY_RATE_$(1)).csv"; \
		if [ ! -f "$$$$src" ]; then \
			echo "  WARNING: $$$$src not found — skipping"; \
			continue; \
		fi; \
		if [ ! -f "$$$$dst" ]; then \
			echo "  Preprocessing $$$${split} target..."; \
			uv run novoboard preprocess \
				--denovo-file "$$$$src" \
				--denovo-output "$$$$dst" \
				--filter-prefix ""; \
		else \
			echo "  $$$$dst already exists — skipping"; \
		fi; \
		if [ -f "$$$$decoy_src" ] && [ ! -f "$$$$decoy_dst" ]; then \
			echo "  Preprocessing $$$${split} decoy (rate=$(DECOY_RATE_$(1)))..."; \
			uv run novoboard preprocess \
				--denovo-file "$$$$decoy_src" \
				--denovo-output "$$$$decoy_dst" \
				--filter-prefix ""; \
		elif [ ! -f "$$$$decoy_src" ]; then \
			echo "  WARNING: $$$$decoy_src not found — skipping"; \
		else \
			echo "  $$$$decoy_dst already exists — skipping"; \
		fi; \
	done
endef
$(foreach origin,$(DATASET_ORIGINS),$(eval $(call preprocess_origin,$(origin))))

preprocess-datasets::
	@echo "Preprocessing complete!"

# Estimate per-PSM FDR and q-values for each origin/split.
# Requires preprocess-datasets to have been run first.
define fdr_origin
# $(1) = origin name
fdr-datasets:: fdr-dataset-$(1)
fdr-dataset-$(1): preprocess-dataset-$(1)
	@echo "=== FDR estimation: $(1) ==="
	@for split in $(DATASET_SPLITS); do \
		target="$(DATASET_DIR)/$(1)/$(DATASET_NOVOBOARD)/$$$${split}.csv"; \
		decoy="$(DATASET_DIR)/$(1)/$(DATASET_NOVOBOARD)/$$$${split}_decoy_$(DECOY_RATE_$(1)).csv"; \
		if [ ! -f "$$$$target" ]; then \
			echo "  WARNING: $$$$target not found — skipping"; \
			continue; \
		fi; \
		if [ ! -f "$$$$decoy" ]; then \
			echo "  WARNING: $$$$decoy not found — skipping"; \
			continue; \
		fi; \
		echo "  Estimating FDR for $$$${split} (decoy rate=$(DECOY_RATE_$(1)))..."; \
		uv run novoboard fdr \
			--target-file "$$$$target" \
			--decoy-files "$$$$decoy" \
			--score-column "ALC (%)" \
			--aa-score-column "local confidence (%)" \
			--tp-metric peptide; \
	done
endef
$(foreach origin,$(DATASET_ORIGINS),$(eval $(call fdr_origin,$(origin))))

fdr-datasets::
	@echo "FDR estimation complete!"
