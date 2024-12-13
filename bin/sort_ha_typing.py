import pandas as pd
import os
import sys

def main(excel_file, fasta_dir):
    # Read the Excel file
    df = pd.read_excel(excel_file, sheet_name="Subtype Predictions")
    
    # Extract relevant columns
    df = df[["Sample", "Subtype Prediction"]]

    # Create a dictionary to hold sequences for each subtype
    subtype_sequences = {
        "H1N1": [],
        "H3N2": [],
        "H5N1": [],
        "VIC": [],
        "YAM": []
    }

    # Group sequences by subtype
    for _, row in df.iterrows():
        sample_id = row["Sample"]
        subtype = row["Subtype Prediction"]

        # Determine the subtype and corresponding file
        if "H1N1" in subtype:
            subtype_file = "H1N1_samples.fasta"
        elif "H3N2" in subtype:
            subtype_file = "H3N2_samples.fasta"
        elif "H5N1" in subtype:
            subtype_file = "H5N1_samples.fasta"
        elif "VIC" in subtype:
            subtype_file = "VIC_samples.fasta"
        elif "YAM" in subtype:
            subtype_file = "YAM_samples.fasta"
        else:
            continue

        # Read the sample FASTA file and append its contents to the corresponding subtype list
        sample_fasta = os.path.join(fasta_dir, f"{sample_id}.fasta")
        if os.path.exists(sample_fasta):
            with open(sample_fasta, 'r') as f:
                sequence = f.read()
                subtype_sequences[subtype].append(sequence)
    
    # Write the sequences to multi-FASTA files
    for subtype, sequences in subtype_sequences.items():
        output_file = f"{subtype}_samples.fasta"
        with open(output_file, 'w') as f:
            for seq in sequences:
                f.write(seq)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python sort_h_typing.py <excel_file> <fasta_dir>")
        sys.exit(1)
    
    excel_file = sys.argv[1]
    fasta_dir = sys.argv[2]
    main(excel_file, fasta_dir)
