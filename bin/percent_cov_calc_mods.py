
#!/usr/bin/env python

# version = '1.0.0'

# Modified from script https://github.com/CDPHE-bioinformatics/CDPHE-influenza/blob/main/scripts/calc_percent_cov.py

# import python modules
import pandas as pd
from datetime import date
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import sys
import argparse
import subprocess

### Segment length dictionary
ref_len_dict = {
    "A_MP": 982,
    "A_NP": 1497,
    "A_NS": 863,
    "A_PA": 2151,
    "A_PB1": 2274,
    "A_PB2": 2280,
    "A_HA_H1": 1704,
    "A_HA_H10": 1686,
    "A_HA_H11": 1698,
    "A_HA_H12": 1695,
    "A_HA_H13": 1701,
    "A_HA_H14": 1707,
    "A_HA_H15": 1713,
    "A_HA_H16": 1698,
    "A_HA_H2": 1689,
    "A_HA_H3": 1704,
    "A_HA_H4": 1695,
    "A_HA_H5": 1707,
    "A_HA_H6": 1704,
    "A_HA_H7": 1713,
    "A_HA_H8": 1701,
    "A_HA_H9": 1683,
    "A_NA_N1": 1413,
    "A_NA_N2": 1410,
    "A_NA_N3": 1410,
    "A_NA_N4": 1413,
    "A_NA_N5": 1422,
    "A_NA_N6": 1413,
    "A_NA_N7": 1416,
    "A_NA_N8": 1413,
    "A_NA_N9": 1413,
    "B_HA": 1758,
    "B_MP": 1139,
    "B_NA": 1408,
    "B_NP": 1683,
    "B_NS": 1034,
    "B_PA": 2181,
    "B_PB1": 2263,
    "B_PB2": 2313,
}


#### FUNCTIONS #####
def getOptions():
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("fasta_files", help="Path to the fasta file")
    parser.add_argument("meta_id", help="Meta ID")
    options = parser.parse_args()
    return options


def get_fasta_file_basename(fasta_file_path):
    basename = fasta_file_path.split("/")[-1]  # strip directories
    print(fasta_file_path)
    return basename


def get_segment_name(fasta_file_path):
    basename = fasta_file_path.split("/")[-1]  # strip directories
    segment_name = basename.split(".")[0]  # remove file extension
    return segment_name


def get_gene_name(fasta_file_path):
    basename = fasta_file_path.split("/")[-1]  # strip directories
    segment_name = basename.split(".")[0]  # remove file extension

    gene_name = segment_name.split("_")  # Split using "_"
    
    if len(gene_name) > 1:
        return gene_name[1]  # Extract gene name
    else:
        print(f"Warning: Could not extract gene name from '{basename}'. Skipping.")
        return None
    
# def get_gene_name(fasta_file_path):
#     basename = fasta_file_path.split("/")[-1]  # strip directories
#     segment_name = basename.split(".")[0]  # remove file extension
#     gene_name = segment_name.split("_")[1]  # extract gene name
#     if len(gene_name) > 1:
#         return gene_name[1]  # extract gene name
#     else:
#         return None  # or raise an error


def get_seq_length(fasta_file_path):
    # read in fasta file
    record = SeqIO.read(fasta_file_path, "fasta")

    # get length of non ambigous bases
    seq = record.seq
    seq_length = seq.count("A") + seq.count("C") + seq.count("G") + seq.count("T")

    return seq_length


def calc_percent_cov(seq_length, ref_len_dict, segment_name):
    if segment_name not in ref_len_dict:
        print(f"Warning: Segment name '{segment_name}' not found in reference dictionary.")
        return None  # Handle the case properly    # calcuate per cov based on expected ref length

    expected_length = ref_len_dict[segment_name]
    percent_coverage = round(((seq_length / expected_length) * 100), 2)

    return percent_coverage


def create_output(meta_id, segment_name, seq_length, percent_coverage, reference_length):
    df = pd.DataFrame()
    df["Sample"] = [meta_id]
    df["segment_name"] = [segment_name]
    df["reference_length"] = [reference_length]
    df["seq_length"] = [seq_length]
    df["percent_coverage"] = [percent_coverage]

    # Construct the output filename header
    output_header = "Sample\tsegment_name\treference_length\tseq_length\tpercent_coverage"

    # Construct the output filename
    output_filename = f'{meta_id}.{segment_name}.perc_cov_results.tsv'

    # Write the dataframe to the output file
    with open(output_filename, "w") as f:
        f.write(output_header + "\n")
        df.to_csv(f, sep="\t", index=False, header=False)


#### MAIN ####
if __name__ == "__main__":
    options = getOptions()
    fasta_file_path = options.fasta_files
    meta_id = options.meta_id

    basename = get_fasta_file_basename(fasta_file_path=fasta_file_path)
    
    segment_name = get_segment_name(fasta_file_path=fasta_file_path)
    gene_name = get_gene_name(fasta_file_path=fasta_file_path)
    
    seq_length = get_seq_length(fasta_file_path=fasta_file_path)
    percent_coverage = calc_percent_cov(seq_length=seq_length, ref_len_dict=ref_len_dict, segment_name=segment_name)
    
    reference_length = ref_len_dict[segment_name]

 # **Skip processing if gene_name is missing**
    if not (".fa" in fasta_file_path or ".fasta" in fasta_file_path):
        print(f"Warning: File '{fasta_file_path}' does not have a valid FASTA extension (.fa or .fasta). Skipping.")
    
    if gene_name is None:
        print("Warning: Could not extract gene name from filename. Skipping this file.")
    
    if seq_length is None:
        print("Warning: Could not calculate seq length. Skipping this file.")
    
    if percent_coverage is None:
        print(f"Warning: Could not calculate percent coverage for segment '{segment_name}'. Skipping.")
    
    if segment_name not in ref_len_dict:
        print(f"Warning: Segment name '{segment_name}' not found in reference dictionary. Skipping.")
        

    create_output(
        meta_id=meta_id,
        segment_name=segment_name,
        seq_length=seq_length,
        percent_coverage=percent_coverage,
        reference_length=reference_length,
        )
    # seq_length = get_seq_length(fasta_file_path=fasta_file_path)
    # percent_coverage = calc_percent_cov(seq_length=seq_length, ref_len_dict=ref_len_dict, segment_name=segment_name)

    # if gene_name is None:
    #     print("Error: Could not extract gene name from filename. Please check the input file format.")
    #     sys.exit(1)

    # if percent_coverage is None:
    #     print("Error: Could not calculate percent coverage. Please check the segment name.")
    #     sys.exit(1)

    # reference_length = ref_len_dict[segment_name]

    # create_output(
    #     meta_id=meta_id,
    #     segment_name=segment_name,
    #     seq_length=seq_length,
    #     percent_coverage=percent_coverage,
    #     reference_length=reference_length,
    # )
