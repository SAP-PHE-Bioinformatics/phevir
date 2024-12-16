#!/usr/bin/env python


import pandas as pd
import argparse
from sqlalchemy import create_engine
from datetime import datetime
import os
#from dotenv import load_dotenv

def process_covid_data(pandepth_path, pangolin_path, nextclade_path, pangocollapse_path, run_id, fasta_dir):
    # Read pandepth file
    # df_pandepth = pd.read_csv(pandepth_path, sep='\t')
    # # drop any row starting with ## as they are comments
    # df_pandepth = df_pandepth[~df_pandepth.apply(lambda row: row.astype(str).str.startswith('##').any(), axis=1)]
    # df_pandepth['Chr'] = df_pandepth['Chr'].astype(str)
    # df_pandepth['Sample'] = df_pandepth['Chr'].astype(str)
    df_pandepth = pd.read_csv(
    pandepth_path,
    sep='\t',
    header=0,
    names=['Sample','Chr', 'Length', 'CoveredSite', 'TotalDepth', 'Coverage(%)', 'MeanDepth'],
    comment='##',  # Skip rows starting with '#'
    engine='python'  # Use Python engine to handle mixed headers and data
)   
    df_pandepth.rename(columns={'Coverage(%)': "Coverage",}, inplace=True)
    df_pandepth['Sample'] = df_pandepth['Sample'].astype(str)
    df_pandepth.head()

    # Read pangolin file
    df_pangolin = pd.read_csv(pangolin_path, sep=',')
    df_pangolin['Sample'] = df_pangolin['taxon'].str.split('/').str[0]

    # Read nextclade file
    df_nextclade = pd.read_csv(nextclade_path, sep='\t')
    df_nextclade['Sample'] = df_nextclade['seqName'].str.split('/').str[0]

    # Read pangocollapse file
    df_pangocollapse = pd.read_csv(pangocollapse_path, sep='\t')
    #sample name lloks like SA2178835/ARTIC/medaka need to drop the /ARTIC/medaka so do follwing
    df_pangocollapse['Sample'] = df_pangocollapse['seqName'].str.split('/').str[0]


    # Merge all dataframes on 'Sample'
    df_merged = pd.merge(df_pandepth, df_pangolin[['Sample', 'lineage']], on='Sample', how='left')
    df_merged = pd.merge(df_merged, df_nextclade[['Sample', 'Nextclade_pango']], on='Sample', how='left')
    df_merged = pd.merge(df_merged, df_pangocollapse[['Sample', 'clade_who', 'VOC_Lineage']], on='Sample', how='left')

    # Add RunID and RunDate
    df_merged['RunID'] = f'{run_id}_COV'
    df_merged['date'] = datetime.now().strftime('%Y-%m-%d')

    df_merged['MeanDepth'] = pd.to_numeric(df_merged['MeanDepth'], errors='coerce')
    df_merged['Coverage'] = pd.to_numeric(df_merged['Coverage'], errors='coerce')

    # Fill NaN values with 0
    df_merged['MeanDepth'] = df_merged['MeanDepth'].fillna(0)
    df_merged['Coverage'] = df_merged['Coverage'].fillna(0)
    df_merged.rename(columns={'Coverage': '%GenomeFraction', 'MeanDepth': 'MeanDepthCov'}, inplace=True)

    print(df_merged.head())


    # Add QC column
    df_merged['GenomicQC'] = df_merged.apply(
        lambda row: 'PASS' if row['MeanDepthCov'] > 100 and row['%GenomeFraction'] > 90 else 'FAIL', axis=1
    )
    # if GenomicQC pass, then set Species_Identification to SARS-CoV-2
    df_merged['Species_Identification'] = df_merged.apply(
        lambda row: 'SARS-CoV-2' if row['GenomicQC'] == 'PASS' else '', axis=1
    )

    # Reorder columns for clarity
    columns_order = ['Sample', 'RunID', 'date', 'MeanDepthCov', '%GenomeFraction', 'GenomicQC', 'clade_who', 'Nextclade_pango'
                     , 'VOC_Lineage', 'lineage']
    df_merged = df_merged[columns_order]

    # Output to CSV
    output_csv_path = f"{run_id}_covid_summary.csv"
    df_merged.to_csv(output_csv_path, index=False)
    print(f"Summary CSV saved to {output_csv_path}")

    def concatenate_fasta_files(fasta_files, output_path, format_type="default", dataset_name=None):
        """
        Concatenates multiple FASTA files into a single file.

        Parameters:
        - fasta_files: List of FASTA file paths to concatenate.
        - output_path: Output path for the concatenated FASTA file.
        - format_type: "default" for standard concatenation, "gisaid" for GISAID-specific formatting.
        - dataset_name: Optional name for the dataset (used in GISAID headers).
        """
        current_year = datetime.now().strftime('%Y')
        with open(output_path, 'w') as outfile:
            for fasta_file in fasta_files:
                with open(fasta_file, 'r') as infile:
                    for line in infile:
                        if line.startswith('>'):  # Header line
                            original_header = line[1:].strip()  # Remove '>' and strip whitespace
                            sample_name = original_header.split('/')[0]
                            if format_type == "gisaid":
                                # Format header for GISAID
                                header = f">hCoV-19/Australia/{sample_name}/{current_year}\n"
                                outfile.write(header)
                            else:
                                # Keep the original header for default
                                outfile.write(line)
                        else:
                            # Write sequence lines as-is
                            outfile.write(line)

        print(f"Concatenated FASTA file saved to {output_path}")

    # Filter samples that pass QC
    passed_samples = df_merged[df_merged['GenomicQC'] == 'PASS']
    passed_samples_path = f"{run_id}_passed_samples.txt"
    passed_samples[['Sample']].to_csv(passed_samples_path, index=False, header=False)
    print(f"Passed samples saved to {passed_samples_path}")


    concat_fasta_path = f"{run_id}_concatenated.fasta"
    gisaid_fasta_path = f"{run_id}_gisaid.fasta"

    # Generate concatenated FASTA with sample names in headers
    concatenate_fasta_files(
    fasta_files=[os.path.join(fasta_dir, f"{sample}.consensus.fasta") for sample in passed_samples['Sample']],
    output_path=concat_fasta_path,
    format_type="default"
)

    # Generate GISAID-formatted FASTA
    concatenate_fasta_files(
    fasta_files=[os.path.join(fasta_dir, f"{sample}.consensus.fasta") for sample in passed_samples['Sample']],
    output_path=gisaid_fasta_path,
    format_type="gisaid",
    dataset_name=run_id
)
    # # Generate concatenated FASTA for passed samples
    # fasta_files = [os.path.join(fasta_dir, f"{sample}.consensus.fasta") for sample in passed_samples['Sample']]
    # concat_fasta_path = f"{run_id}_concatenated.fasta"
    # concatenate_fasta_files(fasta_files, concat_fasta_path)

    # # Generate GISAID formatted FASTA
    # gisaid_fasta_path = f"{run_id}_gisaid.fasta"
    # concatenate_fasta_files(fasta_files, gisaid_fasta_path, format_type="gisaid", dataset_name=run_id)

    # Generate proforma file
    proforma_path = f"{run_id}_proforma.csv"
    with open(proforma_path, 'w') as proforma_file:
        proforma_file.write("Seq_ID,Owner_group,Shared_groups\n")
        for sample in passed_samples['Sample']:
            proforma_file.write(f"{sample},SAP-Owner,SC2-ANZ-Group;SAP-Everyone\n")
    print(f"Proforma saved to {proforma_path}")

    # create summary csv
    summary_csv = f"{run_id}_covid_summary.csv"
    df_merged.to_csv(summary_csv, index=False)

   

    return output_csv_path, passed_samples_path, concat_fasta_path, gisaid_fasta_path, proforma_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process COVID data and output results.')
    parser.add_argument('--pandepth', required=True, help='Path to the pandepth file.')
    parser.add_argument('--pangolin', required=True, help='Path to the pangolin file.')
    parser.add_argument('--nextclade', required=True, help='Path to the nextclade file.')
    parser.add_argument('--pangocollapse', required=True, help='Path to the pangocollapse file.')
    parser.add_argument('--run', required=True, help='RunID must be entered e.g. COVID23003')
    parser.add_argument('--fasta_dir', required=True, help='Path to the directory containing fasta files.')

    args = parser.parse_args()
    process_covid_data(args.pandepth, args.pangolin, args.nextclade, args.pangocollapse, args.run, args.fasta_dir)
