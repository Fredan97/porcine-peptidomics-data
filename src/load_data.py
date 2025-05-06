import os
import pandas as pd
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.process_data.peaks_processor import process_data, preprocess_data_matrix


def load_main_data(data_dir="data", output_dir="outputs", save_files=True):
    """
    Load and process the main experimental data from Pig_day directories.
    
    Parameters:
    -----------
    data_dir : str
        Directory containing the data files and design.csv
    output_dir : str
        Directory where processed data will be saved
    save_files : bool
        Whether to save the processed data to CSV files
        
    Returns:
    --------
    tuple
        (data_matrix, design_matrix, processed_data_matrix)
    """
    if save_files:
        os.makedirs(output_dir, exist_ok=True)
        
    # Process the main data from Pig_day directories
    data_matrix, design_matrix = process_data(
        data_dir, os.path.join(data_dir, "design.csv")
    )
    processed_data_matrix = preprocess_data_matrix(data_matrix)
    
    if save_files:
        processed_data_matrix.to_csv(os.path.join(output_dir, "processed_data_matrix.csv"))
        data_matrix.to_csv(os.path.join(output_dir, "data_matrix.csv"))
        design_matrix.to_csv(os.path.join(output_dir, "design_matrix.csv"))
    
    return data_matrix, design_matrix, processed_data_matrix


def load_rerun_data(data_dir="data", output_dir="outputs", save_files=True):
    """
    Load and process the rerun experimental data from Rerun_files directory.
    
    Parameters:
    -----------
    data_dir : str
        Base directory containing Rerun_files subdirectory and rerun_design.csv
    output_dir : str
        Directory where processed data will be saved
    save_files : bool
        Whether to save the processed data to CSV files
        
    Returns:
    --------
    tuple
        (rerun_data_matrix, rerun_design_matrix, processed_rerun_data_matrix)
    """
    if save_files:
        os.makedirs(output_dir, exist_ok=True)
        
    rerun_dir = os.path.join(data_dir, "Rerun_files")
    rerun_design_file = os.path.join(data_dir, "rerun_design.csv")
    
    # Process rerun data
    rerun_data_matrix, rerun_design_matrix = process_data(rerun_dir, rerun_design_file)
    processed_rerun_data_matrix = preprocess_data_matrix(rerun_data_matrix)
    
    if save_files:
        rerun_data_matrix.to_csv(os.path.join(output_dir, "rerun_data_matrix.csv"))
        processed_rerun_data_matrix.to_csv(os.path.join(output_dir, "processed_rerun_data_matrix.csv"))
        rerun_design_matrix.to_csv(os.path.join(output_dir, "rerun_design_matrix.csv"))
    
    return rerun_data_matrix, rerun_design_matrix, processed_rerun_data_matrix