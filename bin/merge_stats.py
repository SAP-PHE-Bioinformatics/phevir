#! /usr/bin/env python
import pandas as pd
import gzip
import sys
import os

def read_pandepth(pandepth_file, sample, segment):
    """Reads and processes the pandepth file, keeping only relevant columns."""
    with open(pandepth_file, 'r') as f:
        lines = f.readlines()
    
    # Extract header and data
    header = "Sample\tChr\tTotalDepth\tMeanDepth\n"
    relevant_data = []
    for line in lines:
        if line.startswith("#") or line.startswith("##"):
            continue  # Skip comment lines
        parts = line.strip().split()
        relevant_data.append(f"{sample}\t{parts[0]}\t{parts[3]}\t{parts[5]}")

    return header, relevant_data

def read_coverage(coverage_file):
    """Reads the coverage file and excludes the first two columns."""
    coverage_df = pd.read_csv(coverage_file, sep="\t")
    coverage_df = coverage_df.iloc[:, 2:]  # Drop first two columns (Sample and segment_name)
    return coverage_df

def merge_files(pandepth_data, coverage_file, output_file):
    """Merges the processed pandepth data with the coverage file."""
    # Write pandepth data to a temporary DataFrame
    pandepth_df = pd.DataFrame(
        [row.split('\t') for row in pandepth_data],
        columns=["Sample", "Segment", "TotalDepth", "MeanDepth"]
    )

    # Read coverage file and merge
    coverage_df = read_coverage(coverage_file)
    merged_df = pd.concat([pandepth_df, coverage_df], axis=1)
    
    # Save to output
    merged_df.to_csv(output_file, sep="\t", index=False)

def main():
    if len(sys.argv) != 5:
        print("Usage: python merge_stats.py <pandepth_file> <coverage_file> <sample> <segment>")
        sys.exit(1)
    
    pandepth_file = sys.argv[1]
    coverage_file = sys.argv[2]
    sample = sys.argv[3]
    segment = sys.argv[4]
    
    prefix = sample + "_" + segment
    output_file = f"{prefix}_merged.stats"
    
    pandepth_header, pandepth_data = read_pandepth(pandepth_file, sample, segment)
    if not pandepth_data:
        print(f"No matching data found for sample: {sample}, segment: {segment} in {pandepth_file}")
        sys.exit(1)
    
    # Write processed pandepth data to a file
    with open(f"{prefix}_pandepth.stats", "w") as f:
        f.write(pandepth_header + "\n".join(pandepth_data) + "\n")
    
    merge_files(pandepth_data, coverage_file, output_file)
    print(f"Merged file saved as: {output_file}")

if __name__ == "__main__":
    main()