import os
import argparse
import pandas as pd
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.process_data.peaks_processor import process_data, preprocess_data_matrix
from src.load_data import load_main_data, load_rerun_data


def main():
    """
    Process peptidomics data from PEAKS software output and save results to CSV files.
    
    This script processes both main experimental data and rerun data, applying 
    appropriate preprocessing steps for mass spectrometry peptidomics data.
    """
    parser = argparse.ArgumentParser(description="Process peptidomics data from PEAKS software output")
    parser.add_argument("--data-dir", default="data", help="Directory containing the data files")
    parser.add_argument("--output-dir", default="outputs", help="Directory to save processed data")
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f"Processing main experimental data from {args.data_dir}...")
    data_matrix, design_matrix, processed_data_matrix = load_main_data(
        data_dir=args.data_dir, 
        output_dir=args.output_dir
    )
    print(f"Main data processing complete. Processed {data_matrix.shape[1]} samples and {data_matrix.shape[0]} peptides.")
    

    print(f"Processing rerun data from {os.path.join(args.data_dir, 'Rerun_files')}...")
    rerun_data_matrix, rerun_design_matrix, processed_rerun_data_matrix = load_rerun_data(
        data_dir=args.data_dir, 
        output_dir=args.output_dir
    )
    print(f"Rerun data processing complete. Processed {rerun_data_matrix.shape[1]} samples and {rerun_data_matrix.shape[0]} peptides.")

    print(f"All data has been processed and saved to {args.output_dir}/")


if __name__ == "__main__":
    main()
