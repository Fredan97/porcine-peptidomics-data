# Porcine Peptidomics Data

A repository for processing and analyzing mass spectrometry peptidomics data from infected and uninfected porcine wounds.

## Project Overview

This repository contains data and analysis code for mass spectrometry-based peptidomics analysis of porcine wound fluid. The project includes data from the main experiment conducted over three days and a rerun experiment for validation.

### Data Description

- **Main Data**: Located in `data/Pig_day_*` directories, containing PEAKS software exports for each sample
- **Rerun Data**: Located in `data/Rerun_files/`, containing additional validation samples
- **Design Files**: 
  - `data/design.csv`: Contains metadata about main experiment samples
  - `data/rerun_design.csv`: Contains metadata about rerun experiment samples

More information can be found in [/data](/data/README.md)

### Processed Data

All processed data files are stored in the `outputs/` directory:
- `data_matrix.csv`: Raw peptide intensities from main experiment
- `processed_data_matrix.csv`: Log-transformed and preprocessed main data
- `rerun_data_matrix.csv`: Raw peptide intensities from rerun experiment
- `processed_rerun_data_matrix.csv`: Log-transformed and preprocessed rerun data
- `design_matrix.csv` and `rerun_design_matrix.csv`: Sample metadata

## Data Search Parameters

The mass spectrometry data files were searched using PEAKS software with the following parameters:
- Database: Sus scrofa (pig) UniProt database
- Enzyme: No enzyme (peptidomics)
- Precursor mass tolerance: 10 ppm
- Fragment mass tolerance: 0.02 Da
- Variable modifications: Oxidation (M), Acetylation (N-term)

The complete dataset has been deposited to ProteomeXchange with the identifier PXD048892. The blinded re-run dataset is deposited with the identifier PXD055074. 

## Repository Structure

```
data/                   # Raw data files from PEAKS
  design.csv            # Main experiment design file
  rerun_design.csv      # Rerun experiment design file
  Pig_day_1/            # Main experiment day 1 data
  Pig_day_2/            # Main experiment day 2 data 
  Pig_day_3/            # Main experiment day 3 data
  Rerun_files/          # Rerun experiment data
outputs/                # Processed data files
  figures/              # Generated plots and figures

src/                    # Source code
  load_data.py          # Functions to load data matrices
  process_peptides.py   # Script to process raw data
  plot_figures.py       # Script to generate figures
  generate_plots/       # Modules for figure generation
  process_data/         # Modules for data processing
```


## Usage

### Data Processing

To process the raw PEAKS data and generate data matrices:

```bash
python src/process_peptides.py --data-dir data --output-dir outputs
```

Options:
- `--data-dir`: Directory containing the data files (default: "data")
- `--output-dir`: Directory to save processed data (default: "outputs")

### Data Loading API

You can use the data loading API in your scripts to access the processed data:

```python
from src.load_data import load_main_data, load_rerun_data

# Load main experimental data
data_matrix, design_matrix, processed_data_matrix = load_main_data()

# Load rerun data
rerun_data_matrix, rerun_design_matrix, processed_rerun_data_matrix = load_rerun_data()
```

### Generating Figures

To generate figures from the processed data:

```bash
python src/plot_figures.py
```
