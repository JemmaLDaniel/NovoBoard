# NovoBoard

A comprehensive framework for evaluating the false discovery rate and accuracy of de novo peptide sequencing.

## Features

- Calculate fragment ion, amino acid, and peptide accuracies
- Generate decoy spectra for FDR estimation
- Validate FDR estimation against known database matches
- Visualization of FDR validation results
- **Software-agnostic**: Works with any de novo sequencing tool

## Installation

### Using uv (recommended)

```bash
# Clone the repository
git clone https://github.com/BioGeek/NovoBoard.git
cd NovoBoard

# Create virtual environment and install dependencies
uv sync

# Activate the virtual environment
source .venv/bin/activate
```

### Using pip

```bash
pip install -e .
```

## Usage

NovoBoard provides a command-line interface with three main commands:

### 1. Calculate Accuracy

Compare de novo predictions against database search results:

```bash
novoboard accuracy \
    --db-file db_results.csv \
    --denovo-file denovo_predictions.csv \
    --spectrum-file spectra.mgf \
    --score-column "ALC (%)" \
    --aa-score-column "local confidence (%)"
```

**Options:**
- `--db-file`: Database search results CSV (required)
- `--denovo-file`: De novo sequencing results CSV (required)
- `--spectrum-file`: MGF spectrum file (required)
- `--score-column`: Column name for peptide score (default: "ALC (%)")
- `--aa-score-column`: Column name for AA-level scores (default: "local confidence (%)")

**Output:**
- Files are saved in the same directory as the de novo file
- `{denovo_stem}_accuracy.csv`: Per-peptide accuracy metrics
- `{denovo_stem}_denovo_only.csv`: Peptides found only in de novo results
- `{denovo_stem}_scan2fea.csv`: Scan to feature mapping
- `{denovo_stem}_multifea.csv`: Multi-feature entries

Example: `results/denovo.csv` → `results/denovo_accuracy.csv`, etc.

### 2. Generate Decoy Spectra

Create decoy MGF files for FDR estimation:

```bash
novoboard decoy \
    --spectrum-file spectra.mgf \
    --sampling-strategy random \
    --sampling-rate 0.5 \
    --seed 99
```

**Output:**
- Files are saved in the same directory as the input
- Naming convention: `{input_stem}_decoy_{rate}.mgf`
- Example: `spectra.mgf` → `spectra_decoy_0.50.mgf`
- For `permutation` strategy: `spectra_permutation.mgf`
- For `500Da` strategy: `spectra_500Da.mgf`

**Options:**
- `--spectrum-file`: Input MGF file(s) (required, accepts multiple)
- `--sampling-strategy`: Strategy for peak sampling (default: random)
  - `random`: Random peak sampling
  - `intensity`: Sample by peak intensity
  - `intensity_mass`: Intensity + peptide mass
  - `permutation`: Shuffle intensities
  - `500Da`: Remove peaks < 500 Da
  - `distance`: Remove by AA mass distances
- `--sampling-rate`: Fraction of peaks to keep (default: 0.5)
- `--seed`: Random seed for reproducibility (default: 99)

### 3. Validate FDR

Validate FDR estimation using target-decoy approach:

```bash
novoboard fdr \
    --target-file target_denovo.csv \
    --decoy-files decoy1.csv decoy2.csv decoy3.csv \
    --db-file db_results.csv \
    --spectrum-file spectra.mgf \
    --output-file fdr_validation.png \
    --ion-threshold 0.90
```

**Options:**
- `--target-file`: Target de novo results CSV (required)
- `--decoy-files`: Decoy de novo results CSV(s) (required, accepts multiple)
- `--db-file`: Database search results CSV (required)
- `--spectrum-file`: MGF spectrum file (required)
- `--output-file`: Output plot path (default: fdr_validation.png)
- `--score-column`: Column name for peptide score (default: "ALC (%)")
- `--aa-score-column`: Column name for AA-level scores (default: "local confidence (%)")
- `--ion-threshold`: Ion matching threshold (default: 0.90)
- `--labels`: Custom labels for each decoy file (in same order as `--decoy-files`)
- `--fdr-max`: Maximum FDR value for plot axes, between 0 and 1 (default: 0.05)
- `--dpi`: Plot resolution in dots per inch (default: 150, use 300 for print quality)
- `--no-monotonic`: Disable monotonic filtering to show all FDR data points

**Labels:**

If `--labels` is not provided, labels are automatically extracted from filenames. For example, `helaqc_decoy_0.30_results.csv` becomes `30%`.

```bash
# With explicit labels (must match order of --decoy-files)
novoboard fdr \
    --target-file target.csv \
    --decoy-files decoy_10.csv decoy_20.csv decoy_30.csv \
    --labels "10%" "20%" "30%" \
    --db-file db.csv \
    --spectrum-file spectra.mgf

# Without labels (auto-extracted from filenames)
novoboard fdr \
    --target-file target.csv \
    --decoy-files decoy_0.10.csv decoy_0.20.csv decoy_0.30.csv \
    --db-file db.csv \
    --spectrum-file spectra.mgf
```

### 4. Download Example Data

Download the ABRF example dataset:

```bash
novoboard download --output-dir data
```

### 5. Preprocess (Convert InstaNovo Output)

Convert InstaNovo predictions to NovoBoard format:

```bash
# Convert de novo predictions
novoboard preprocess \
    --denovo-file data/instanovo_preds.csv \
    --denovo-output results/denovo.csv

# Convert both de novo predictions and labeled MGF (database annotations)
novoboard preprocess \
    --denovo-file data/instanovo_preds.csv \
    --denovo-output results/denovo.csv \
    --db-mgf-file data/labeled_spectra.mgf \
    --db-output results/db_results.csv
```

**Options:**
- `--denovo-file`: Path to InstaNovo predictions CSV
- `--denovo-output`: Path for converted de novo results CSV
- `--db-mgf-file`: Path to labelled MGF file (for database annotations)
- `--db-output`: Path for database results CSV

**Notes:**
- Converts UNIMOD notation to NovoBoard format (e.g., `C[UNIMOD:4]` → `C(+57.02)`)
- Filters peptides with unsupported modifications
- Converts InstaNovo log probabilities to 0-100 scale scores

## Example Workflow

```bash
# 1. Download example data
novoboard download --output-dir data

# 2. Generate decoy spectra
novoboard decoy --spectrum-file data/spectra.mgf

# 3. Run de novo sequencing on target and decoy spectra
# (using your preferred tool: PEAKS, Casanovo, InstaNovo, etc.)

# 4. (If using InstaNovo) Preprocess predictions to NovoBoard format
novoboard preprocess \
    --denovo-file instanovo_output/predictions.csv \
    --denovo-output results/denovo.csv \
    --db-mgf-file instanovo_input/spectra.mgf \
    --db-output results/db_output.csv

# 5. Calculate accuracy
novoboard accuracy \
    --db-file data/db_results.csv \
    --denovo-file results/denovo.csv \
    --spectrum-file data/spectra.mgf

# 6. Validate FDR
novoboard fdr \
    --target-file results/target.csv \
    --decoy-files results/decoy*.csv \
    --db-file data/db_results.csv \
    --spectrum-file data/spectra.mgf \
    --output-file figures/fdr.png
```

## Jupyter Notebook

The original notebook `aa.fdr_github.ipynb` is still available for interactive analysis:

```bash
source .venv/bin/activate
jupyter notebook aa.fdr_github.ipynb
```

## Development

### Running Tests

```bash
source .venv/bin/activate
pytest tests/ -v
```

### Project Structure

```
novoboard/
├── src/novoboard/
│   ├── __init__.py
│   ├── cli.py           # Command-line interface
│   ├── config.py        # Vocabulary and mass definitions
│   ├── accuracy.py      # Accuracy calculation (WorkerTest)
│   ├── decoy.py         # Decoy MGF generation
│   ├── fdr.py           # FDR calculation and validation
│   ├── mgf.py           # MGF file parser
│   └── plotting.py      # Visualization functions
├── tests/               # Unit tests
├── data/                # Data directory (not in git)
├── aa.fdr_github.ipynb  # Original Jupyter notebook
├── config.py            # Compatibility shim for notebook
├── download_data.py     # Data download script
└── pyproject.toml       # Project configuration
```

## Input File Formats

### De novo Results CSV

Required columns:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `Source File` | `str` | Source spectrum file name (`.mgf` suffix automatically stripped) | `"sample"` or `"sample.mgf"` |
| `Scan` | `int` | Scan number | `1234` |
| `Peptide` | `str` | Peptide sequence with modifications | `"PEPTC(+57.02)DE"` |
| Score column | `float` | Peptide-level confidence score (configurable, default: `ALC (%)`) | `85.5` |
| AA score column | `str` | Comma-separated per-residue confidence scores (configurable, default: `local confidence (%)`) | `"90,85,88,92,87,91,89"` |

**Note:** The `Source File` value is used to match spectra between files. The `.mgf` suffix is automatically stripped when building internal feature IDs (e.g., `sample.mgf` becomes `sample||1234`).

Extra columns will not be included in final results.

### Supported Modifications

NovoBoard supports a **limited set of post-translational modifications**. Peptides containing unsupported modifications will be **silently skipped** during accuracy calculations.

| Residue | Modification | Input Format | Internal Representation | Mass Delta |
|---------|--------------|--------------|-------------------------|------------|
| C | Carbamidomethylation | `C(+57.02)` | `C(Carbamidomethylation)` | +57.02 Da |
| M | Oxidation | `M(+15.99)` | `M(Oxidation)` | +15.99 Da |
| N | Deamidation | `N(+0.98)` | `N(Deamidation)` | +0.98 Da |
| Q | Deamidation | `Q(+0.98)` | `Q(Deamidation)` | +0.98 Da |
| S | Phosphorylation | `S(+79.97)` | `S(Phosphorylation)` | +79.97 Da |
| T | Phosphorylation | `T(+79.97)` | `T(Phosphorylation)` | +79.97 Da |
| Y | Phosphorylation | `Y(+79.97)` | `Y(Phosphorylation)` | +79.97 Da |

**Important limitations:**
- **N-terminal modifications** (e.g., acetylation, carbamylation) are **not supported**
- **Other PTMs** (e.g., methylation, ubiquitination) are **not supported**
- Modification masses must match **exactly** (e.g., `+57.02`, not `+57.021`)
- UNIMOD notation (e.g., `C[UNIMOD:4]`) is **not directly supported** — requires preprocessing

#### Using Other De Novo Tools

If your de novo tool uses a different modification format, you'll need to convert it:

| Tool | Native Format | Conversion Needed |
|------|---------------|-------------------|
| PEAKS | `C(+57.02)` | ✅ Native support |
| Casanovo | `C[UNIMOD:4]` | Convert to `C(+57.02)` |
| InstaNovo | `C[UNIMOD:4]` | Convert to `C(+57.02)` |

Example conversion from UNIMOD to NovoBoard format:

| UNIMOD | NovoBoard |
|--------|-----------|
| `C[UNIMOD:4]` | `C(+57.02)` |
| `M[UNIMOD:35]` | `M(+15.99)` |
| `N[UNIMOD:7]` | `N(+0.98)` |
| `Q[UNIMOD:7]` | `Q(+0.98)` |
| `S[UNIMOD:21]` | `S(+79.97)` |
| `T[UNIMOD:21]` | `T(+79.97)` |
| `Y[UNIMOD:21]` | `Y(+79.97)` |

### Database Search Results CSV

Required columns:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `Source File` | `str` | Source spectrum file name (`.mgf` suffix automatically stripped) | `"sample"` or `"sample.mgf"` |
| `Scan` | `int` | Scan number | `1234` |
| `Peptide` | `str` | Peptide sequence with modifications | `"PEPTC(+57.02)DE"` |

**Note:** The `Source File` value must match the de novo results file for spectrum matching. The `.mgf` suffix is automatically stripped.

Extra columns will not be included in final results.

### MGF Spectrum File

Standard MGF format with:
- `BEGIN IONS` / `END IONS` markers
- `TITLE`, `PEPMASS`, `CHARGE`, `SCANS` headers
- Peak list as `m/z intensity` pairs

Extra headers will not be included in final results.

## Citation

If you use NovoBoard in your research, please cite:

> "NovoBoard: a comprehensive framework for evaluating the false discovery rate and accuracy of de novo peptide sequencing"
> https://doi.org/10.1101/2024.04.16.589668
