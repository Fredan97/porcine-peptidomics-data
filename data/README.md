# Data Directory

This directory contains the raw mass spectrometry data files and metadata for the porcine peptidomics project studying wound fluid from infected and uninfected porcine wounds.

## Directory Structure

```
data/
├── design.csv                                    # Main experiment design file
├── rerun_design.csv                             # Rerun experiment design file  
├── Uniprot_Pig_SwissProt_Altered_FIBA&FIBB_2023_05_11.fasta  # Protein database
├── Pig_day_1/                                   # Day 1 sample data
│   ├── 230524_Fredrik_Forsberg_DB_Search_Sample_1_3/
│   ├── 230524_Fredrik_Forsberg_DB_Search_Sample_2_59/
│   └── ... (52 sample directories)
├── Pig_day_2/                                   # Day 2 sample data
│   └── ... (sample directories)
├── Pig_day_3/                                   # Day 3 sample data
│   └── ... (sample directories)
└── Rerun_files/                                 # Rerun experiment data
    ├── Sample_A/
    ├── Sample_B/
    └── ... (12 sample directories)
```

## File Descriptions

### Design Files

#### `design.csv`
Contains metadata for the main experiment samples (114 samples total).

**Columns:**
- `id`: Sample identifier (e.g., "Sample 1 Day 1")
- `group`: Treatment group
  - `S. aureus` (38 samples): Staphylococcus aureus infection
  - `P. aeruginosa` (33 samples): Pseudomonas aeruginosa infection
  - `Ctrl` (26 samples): Control/uninfected samples
  - `Double infection` (8 samples): Dual bacterial infection
  - `Accidental double infection` (8 samples): Unintended dual infection
- `day`: Time point (Day 1, Day 2, or Day 3)
- `pig`: Pig identifier
- `wound`: Wound identifier
- `filename_pride`: Filename for PRIDE repository submission

#### `rerun_design.csv`
Contains metadata for the rerun experiment samples (12 samples total).

**Columns:**
- `blinded_label`: Blinded sample identifier (Sample A-L)
- `sample_name`: Original sample name from main experiment
- `filename_pride`: Filename for PRIDE repository submission

### Protein Database

#### `Uniprot_Pig_SwissProt_Altered_FIBA&FIBB_2023_05_11.fasta`

### Sample Data Directories

#### `Pig_day_1/`, `Pig_day_2/`, `Pig_day_3/`
Each contains individual sample directories for the different days with PEAKS software output files. 

**Directory Naming Convention:**
All sample directories follow the consistent pattern: `Sample_X/` where X is the sample number.

**Sample Directory Structure:**
Each sample directory (regardless of naming convention) contains:
- `DB search psm.csv`: Database search peptide-spectrum matches
- `de novo only peptides.csv`: De novo sequenced peptides
- `peptide.csv`: Identified peptides with quantification
- `protein-peptides.csv`: Protein-peptide mapping
- `proteins.csv`: Protein identifications and quantification

#### `Rerun_files/`
Contains 12 sample directories (`Sample_A/` through `Sample_L/`) with the same file structure as the main experiment samples.

## Data Processing Notes

- **Search Engine**: PEAKS software
- **Database**: Sus scrofa UniProt SwissProt
- **Enzyme**: No enzyme (peptidomics analysis)
- **Precursor mass tolerance**: 10 ppm
- **Fragment mass tolerance**: 0.02 Da
- **Variable modifications**: Oxidation (M), Acetylation (N-term)

## Data Availability

- **Main Dataset**: ProteomeXchange identifier PXD048892
- **Rerun Dataset**: ProteomeXchange identifier PXD055074

## Usage

The raw data files in this directory are processed using the scripts in the `src/` directory to generate the processed data matrices in the `outputs/` directory. See the main project README for details on data processing workflow.