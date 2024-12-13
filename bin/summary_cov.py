import pandas as pd
import argparse
from sqlalchemy import create_engine
from datetime import datetime
import os
from dotenv import load_dotenv

def import_to_sql(df, table_name):
    load_dotenv()
    sql_connection_string = os.getenv("SQLALCHEMY_DATABASE_URL")
    engine = create_engine("mssql+pyodbc://?driver=ODBC+Driver+17+for+SQL+Server")
    try:
        df.to_sql(table_name, engine, if_exists='append', index=False)
        print("Data has been imported successfully.")
    except Exception as error:
        print("Error connecting to SQL Database, check script and if database online. Results not Transferred!", error)

def concatenate_fasta_files(input_files, output_file, format_type="default", dataset_name=None):
    with open(output_file, 'w') as out_f:
        for fasta_file in input_files:
            with open(fasta_file, 'r') as in_f:
                for line in in_f:
                    if line.startswith('>'):
                        # Modify header based on format type
                        if format_type == "gisaid":
                            new_header = f">hCoV-19/Australia/{dataset_name}/2024"
                        else:
                            new_header = f">{os.path.basename(fasta_file).replace('.consensus.fasta', '')}"
                        out_f.write(new_header + '\n')
                    else:
                        out_f.write(line)
    print(f"Concatenated FASTA saved to {output_file}")

def process_covid_data(pandepth_path, pangolin_path, nextclade_path, pangocollapse_path, run_id, fasta_dir):
    # Read pandepth file
    df_pandepth = pd.read_csv(pandepth_path, sep='\t', names=['Sample', 'MeanDepth', 'Coverage'])
    df_pandepth['Sample'] = df_pandepth['Sample'].astype(str)

    # Read pangolin file
    df_pangolin = pd.read_csv(pangolin_path, sep='\t')
    df_pangolin.rename(columns={'lineage': 'Pangolin_Clade'}, inplace=True)
    df_pangolin['Sample'] = df_pangolin['taxon'].str.extract(r'(\d+)')[0]

    # Read nextclade file
    df_nextclade = pd.read_csv(nextclade_path, sep='\t')
    df_nextclade.rename(columns={'clade': 'Raw_Nextclade_Clade'}, inplace=True)
    df_nextclade['Sample'] = df_nextclade['seqName'].str.extract(r'(\d+)')[0]

    # Read pangocollapse file
    df_pangocollapse = pd.read_csv(pangocollapse_path, sep='\t')
    df_pangocollapse.rename(columns={'collapsed_lineage': 'Collapsed_Clade'}, inplace=True)
    df_pangocollapse['Sample'] = df_pangocollapse['Sample'].astype(str)

    # Merge all dataframes on 'Sample'
    df_merged = pd.merge(df_pandepth, df_pangolin[['Sample', 'Pangolin_Clade']], on='Sample', how='left')
    df_merged = pd.merge(df_merged, df_nextclade[['Sample', 'Raw_Nextclade_Clade']], on='Sample', how='left')
    df_merged = pd.merge(df_merged, df_pangocollapse[['Sample', 'Collapsed_Clade']], on='Sample', how='left')

    # Add RunID and RunDate
    df_merged['RunID'] = run_id
    df_merged['RunDate'] = datetime.now().strftime('%Y-%m-%d')

    # Add QC column
    df_merged['Genomic_QC'] = df_merged.apply(
        lambda row: 'PASS' if row['MeanDepth'] > 100 and row['Coverage'] > 90 else 'FAIL', axis=1
    )

    # Reorder columns for clarity
    columns_order = ['Sample', 'RunID', 'RunDate', 'MeanDepth', 'Coverage', 'Genomic_QC', 
                     'Raw_Nextclade_Clade', 'Collapsed_Clade', 'Pangolin_Clade']
    df_merged = df_merged[columns_order]

    # Output to CSV
    output_csv_path = f"{run_id}_covid_summary.csv"
    df_merged.to_csv(output_csv_path, index=False)
    print(f"Summary CSV saved to {output_csv_path}")

    # Filter samples that pass QC
    passed_samples = df_merged[df_merged['Genomic_QC'] == 'PASS']
    passed_samples_path = f"{run_id}_passed_samples.txt"
    passed_samples[['Sample']].to_csv(passed_samples_path, index=False, header=False)
    print(f"Passed samples saved to {passed_samples_path}")

    # Generate concatenated FASTA for passed samples
    fasta_files = [os.path.join(fasta_dir, f"{sample}.consensus.fasta") for sample in passed_samples['Sample']]
    concat_fasta_path = f"{run_id}_concatenated.fasta"
    concatenate_fasta_files(fasta_files, concat_fasta_path)

    # Generate GISAID formatted FASTA
    gisaid_fasta_path = f"{run_id}_gisaid.fasta"
    concatenate_fasta_files(fasta_files, gisaid_fasta_path, format_type="gisaid", dataset_name=run_id)

    # Generate proforma file
    proforma_path = f"{run_id}_proforma.csv"
    with open(proforma_path, 'w') as proforma_file:
        proforma_file.write("Seq_ID,Owner_group,Shared_groups\n")
        for sample in passed_samples['Sample']:
            proforma_file.write(f"{sample},SAP-Owner,SC2-ANZ-Group;SAP-Everyone\n")
    print(f"Proforma saved to {proforma_path}")

    # Import to SQL
    import_to_sql(df_merged, 'covid_results')

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
