import os
import glob
import pandas as pd
import numpy as np
from typing import Dict, Tuple
import re


def load_peptide_file(filepath: str) -> pd.DataFrame:
    """
    Load a single peptide.csv file exported from PEAKS.
    """
    try:
        df = pd.read_csv(filepath)

        folder_path = os.path.dirname(filepath)

        if "Rerun_files" in filepath:

            sample_folder = os.path.basename(os.path.dirname(filepath))
            sample_id = sample_folder.replace("_", " ")
            print(f"Rerun file: {filepath} matches {sample_id}")
        else:
            day_match = re.search(r"Pig_day_(\d+)", filepath, re.IGNORECASE)
            day = f"Day {day_match.group(1)}" if day_match else None

            sample_match = re.search(r"Sample_(\d+)(?:_\d+)?", filepath)
            if sample_match:
                sample_number = sample_match.group(1)
                sample_id = f"Sample {sample_number} {day}"
                print(f"{filepath} matches sample {sample_id}")
            else:
                dir_name = os.path.basename(os.path.dirname(filepath))
                sample_id = f"{dir_name} {day}"
                print(
                    f"No sample number found in {filepath}, using directory name: {sample_id}"
                )

        df["Sample_ID"] = sample_id
        df["Folder"] = folder_path

        if "Source File" in df.columns and len(df) > 0:
            source_file = df["Source File"].iloc[0]
            df["Source_File"] = source_file

        return df

    except Exception as e:
        print(f"Error loading file {filepath}: {e}")
        return pd.DataFrame()


def load_peptide_files(
    data_dir: str,
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, str], Dict[str, str]]:
    """
    Load all peptide.csv files from a directory.
    """
    if "Rerun_files" in data_dir:
        pattern = os.path.join(data_dir, "Sample_*", "**", "peptide.csv")
    else:
        pattern = os.path.join(data_dir, "Pig_day_*", "**", "peptide.csv")

    peptide_files = glob.glob(pattern, recursive=True)

    if not peptide_files:
        print(f"Warning: No peptide files found using pattern: {pattern}")
    else:
        print(f"Found {len(peptide_files)} peptide files")

    results = {}
    sample_folders = {}
    sample_sources = {}

    for filepath in peptide_files:
        df = load_peptide_file(filepath)
        if not df.empty:
            sample_id = df["Sample_ID"].iloc[0]
            results[sample_id] = df
            sample_folders[sample_id] = df["Folder"].iloc[0]
            if "Source_File" in df.columns:
                sample_sources[sample_id] = df["Source_File"].iloc[0]

    # Return the data and the metadata information
    return results, sample_folders, sample_sources


def create_data_matrix(
    peptide_data: Dict[str, pd.DataFrame], intensity_column_default: str = "Area"
) -> pd.DataFrame:
    """
    Create a data matrix from peptide data with peptides as rows and samples as columns.
    """
    all_peptides = set()
    for df in peptide_data.values():
        all_peptides.update(df["Peptide"].unique())

    data_matrix = pd.DataFrame(index=sorted(all_peptides))

    for sample_id, df in peptide_data.items():
        peptide_intensities = {}
        intensity_column = [
            col for col in df.columns if intensity_column_default in col
        ][0]

        if not intensity_column:
            print(
                f"Warning: Could not find any intensity column for {sample_id}. Using '{intensity_column_default}' as fallback."
            )
            intensity_column = intensity_column_default

        for _, row in df.iterrows():
            peptide = row["Peptide"]
            current_intensity = peptide_intensities.get(peptide, 0)
            intensity = row.get(intensity_column, 0)
            peptide_intensities[peptide] = max(current_intensity, intensity)

        sample_data = pd.Series(peptide_intensities, name=sample_id)
        data_matrix = data_matrix.join(sample_data, how="left")

    return data_matrix


def load_design_matrix(filepath: str) -> pd.DataFrame:
    """
    Load the experimental design matrix.
    """
    try:
        design = pd.read_csv(filepath)

        if "rerun_design" in filepath:
            design = design.set_index("sample_name")
            return design
        elif "id" in design.columns:
            design = design.set_index("id")

        return design
    except Exception as e:
        print(f"Error loading design matrix {filepath}: {e}")
        return pd.DataFrame()


def merge_data_with_design(
    data_matrix: pd.DataFrame, design_matrix: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Merge the data matrix with the design matrix based on sample IDs.
    """
    # Get the samples in both matrices
    data_samples = set(data_matrix.columns)
    design_samples = set(design_matrix.index)
    common_samples = data_samples.intersection(design_samples)

    if not common_samples:
        print("No matching samples found between data and design matrices")
        return data_matrix, design_matrix

    filtered_data = data_matrix[sorted(common_samples)]
    filtered_design = design_matrix.loc[sorted(common_samples)]

    return filtered_data, filtered_design


def process_data(data_dir: str, design_file: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process all peptide data and merge with experimental design.

    """
    peptide_data, sample_folders, sample_sources = load_peptide_files(data_dir)
    data_matrix = create_data_matrix(peptide_data)
    design_matrix = load_design_matrix(design_file)

    design_matrix = design_matrix.copy()

    if "Rerun_files" in data_dir:
        print("Processing rerun files - mapping blinded labels to sample names")
        design_df = pd.read_csv(design_file)
        blinded_map = dict(zip(design_df["blinded_label"], design_df["sample_name"]))

        new_columns = {}
        folder_map = {}
        source_map = {}

        for col in data_matrix.columns:
            if col in blinded_map:
                new_name = blinded_map[col]
                new_columns[col] = new_name

                if col in sample_folders:
                    folder_map[new_name] = sample_folders[col]
                if col in sample_sources:
                    source_map[new_name] = sample_sources[col]
                print(f"Mapping {col} to {new_name}")
            else:
                new_columns[col] = col
                print(f"No mapping found for {col}")

        data_matrix = data_matrix.rename(columns=new_columns)
        sample_folders = {new_columns.get(k, k): v for k, v in sample_folders.items()}
        sample_folders.update(folder_map)
        sample_sources = {new_columns.get(k, k): v for k, v in sample_sources.items()}
        sample_sources.update(source_map)

    merged_data, filtered_design = merge_data_with_design(data_matrix, design_matrix)

    for sample_id in filtered_design.index:
        if sample_id in sample_folders:
            filtered_design.loc[sample_id, "folder"] = sample_folders[sample_id]
        if sample_id in sample_sources:
            filtered_design.loc[sample_id, "source_file"] = sample_sources[sample_id]

    merged_data = merged_data.replace(0, np.nan)
    return merged_data, filtered_design


def preprocess_data_matrix(
    data_matrix: pd.DataFrame, log_transform: bool = True
) -> pd.DataFrame:
    """
    Preprocess the data matrix with common operations for mass spec data.

    Include normalization.
    """
    matrix = data_matrix.copy()

    if log_transform:
        matrix = matrix.replace(0, np.nan)
        matrix = matrix.apply(lambda x: np.log2(x) if x.name != "Sample_ID" else x)

    return matrix
