#!/usr/bin/env python

import pandas as pd
import pathlib

def parse_kraken_report(file_path):
    """Parse a Kraken2 report file and return a DataFrame with species-level matches."""
    kraken = pathlib.Path(file_path)
    df = pd.read_csv(kraken, sep="\t", header=None, names=['percentage', 'frag1', 'frag2', 'code', 'taxon', 'name'])

    # Convert percentages to float and sort by percentage in descending order
    df['percentage'] = df['percentage'].apply(lambda x: float(x.strip('%')) if isinstance(x, str) else float(x))
    df = df.sort_values(by='percentage', ascending=False)

    # Filter for species-level codes ('S')
    df = df[df['code'].isin(['S'])]

    # Reset the index for easier access
    df = df.reset_index(drop=True)
    return df

def extract_top_matches(df, kraken):
    """Extract the top three species matches from the DataFrame."""
    tempdf = pd.DataFrame()
    d = {}
    t = len(df)

    for i in range(min(t, 3)):
        d.update({
            'Isolate': kraken.parts[-1],  # Use the filename as the isolate identifier
            f"#{i+1} Match": df.loc[i, 'name'].strip(),
            f"%{i+1}": df.loc[i, 'percentage']
        })

    tempdf = pd.DataFrame(data=d, index=[0])
    return tempdf

def process_kraken_reports(input_files, output_file):
    """Process multiple Kraken2 report files and save the top matches to a TSV file."""
    kfiles = input_files.split()
    id_table = pd.DataFrame()

    for k in kfiles:
        df = parse_kraken_report(k)
        tempdf = extract_top_matches(df, pathlib.Path(k))

        if id_table.empty:
            id_table = tempdf
        else:
            id_table = pd.concat([id_table, tempdf], ignore_index=True)

    # Reorder the columns
    cols_list = ['Isolate', '#1 Match', '%1', '#2 Match', '%2', '#3 Match', '%3']
    id_table = id_table.reindex(cols_list, axis='columns')

    # Save to the output TSV file
    id_table.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process Kraken2 reports and extract top species matches.")
    parser.add_argument("-i", "--input", required=True, help="Input Kraken2 report files (space-separated)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file to save the top matches")

    args = parser.parse_args()

    process_kraken_reports(args.input, args.output)
