#!/usr/bin/env python
import pandas as pd
import argparse
from sqlalchemy import create_engine
from datetime import datetime
import io, os, warnings
from dotenv import load_dotenv


warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

def import_to_sql(df, table_name):
    load_dotenv('/scratch/pheed/config_files/')
    sql_connection_string = os.getenv("SQLALCHEMY_DATABASE_URL")
    print(sql_connection_string)
    engine = create_engine("mssql+pyodbc://sa:WGSmaster22!@frgeneseqgpu/PHEED?driver=ODBC+Driver+17+for+SQL+Server")
    try:
        df.to_sql(table_name, engine, if_exists='append', index=False)
        print("Data has been imported successfully.")
    except Exception as error:
        print("Error connecting to SQL Database, check script and if database online. Results not Transferred!", error)


def aggregate_mutations(mutations_file):
    with open(mutations_file, 'r') as file:
        lines = file.readlines()

    data = []
    for line in lines:
        if not line.startswith("Flusurver result") and not line.startswith("Reference"):
            data.append(line.strip().split("\t"))

    columns = [
        'Reference', 'Query', '% AA identity', '% length coverage', 'HA vacc. eff.', 'total # mutations', 'Mutation', 
        'Interestlevel (0-3)', 'Subtype marker', '# global occ', '# country occ', 'Prev. reported effect', 
        'Multiple position effect', 'Structural interaction(s)', 'Effect on glycosylation'
    ]

    # Create a DataFrame
    df_mutations = pd.DataFrame(data, columns=columns)

    # Filter and clean data
    df_mutations = df_mutations[df_mutations['Mutation'] != 'no mutations']
    df_mutations['Interestlevel (0-3)'] = df_mutations['Interestlevel (0-3)'].astype(str)
    df_mutations['Sample'] = df_mutations['Query'].str.extract(r'(\d+)')[0]
    df_mutations['Segment'] = df_mutations['Reference'].str.split('_').str[0]

    # Add Segment and Mutation Info as Key-Value
    df_mutations['Mutation_Info'] = df_mutations.apply(
        lambda row: f"{row['Segment']}: {row['Mutation']}", axis=1
    )

    # Group mutations by Sample, Interestlevel, and combine into key-value lists
    df_grouped = df_mutations.groupby(['Sample', 'Interestlevel (0-3)'])['Mutation_Info'].apply(
        lambda x: ', '.join(x)
    ).reset_index()

    # Pivot table to create columns for each interest level
    df_pivot = df_grouped.pivot(index='Sample', columns='Interestlevel (0-3)', values='Mutation_Info').reset_index()

    # Rename columns to include interest level explicitly
    df_pivot.columns = ['Sample'] + [f"Level_{col}_Mutations" for col in df_pivot.columns[1:]]

    # Reorder columns to start with Level_3_Mutations
    level_order = ['Level_3_Mutations', 'Level_2_Mutations', 'Level_1_Mutations', 'Level_0_Mutations']
    existing_levels = [col for col in level_order if col in df_pivot.columns]
    df_pivot = df_pivot[['Sample'] + existing_levels]

    # Replace NaN with '-' for missing levels
    df_pivot = df_pivot.fillna('-')

    return df_pivot


def collapse_depth_columns(df):
    # Identify HA and NA columns for MeanDepth and percent_coverage
    ha_depth_columns = [col for col in df.columns if col.startswith("Depth_HA_")]
    na_depth_columns = [col for col in df.columns if col.startswith("Depth_NA_")]
    ha_coverage_columns = [col for col in df.columns if col.startswith("percent_coverage_HA_")]
    na_coverage_columns = [col for col in df.columns if col.startswith("percent_coverage_NA_")]

    # Collapse HA and NA depth columns
    if ha_depth_columns:
        df['Depth_HA'] = df[ha_depth_columns].mean(axis=1)
    if na_depth_columns:
        df['Depth_NA'] = df[na_depth_columns].mean(axis=1)

    # Collapse HA and NA percent_coverage columns
    if ha_coverage_columns:
        df['percent_coverage_HA'] = df[ha_coverage_columns].mean(axis=1)
    if na_coverage_columns:
        df['percent_coverage_NA'] = df[na_coverage_columns].mean(axis=1)

    # Drop original subtype-specific columns
    df = df.drop(columns=ha_depth_columns + na_depth_columns + ha_coverage_columns + na_coverage_columns, errors='ignore')

    return df

def main(tsv_file_path, excel_file_path, qc_file_path, kraken2_file_path, drug_sensitivity_path, typing_report_path, irma_qc_path, run_id, mut):
    # Read QC file
    df_qc = pd.read_csv(qc_file_path, sep='\t', names=[
        'Sample', 'Segment', 'TotalDepth', 'MeanDepth', 'reference_length', 'seq_length', 'percent_coverage'
    ])
    df_qc['Sample'] = df_qc['Sample'].astype(str)
    df_qc['Segment'] = df_qc['Segment'].str.extract(r'_(\w+)$')[0]
    df_qc['MeanDepth'] = pd.to_numeric(df_qc['MeanDepth'], errors='coerce')
    df_qc['percent_coverage'] = pd.to_numeric(df_qc['percent_coverage'], errors='coerce')
    df_qc_cleaned = df_qc.dropna(subset=['Sample', 'Segment', 'MeanDepth'])
    df_qc_cleaned = df_qc_cleaned.groupby(['Sample', 'Segment'], as_index=False).agg({'MeanDepth': 'mean', 'percent_coverage': 'mean'})
    df_qc_pivot = df_qc_cleaned.pivot(index='Sample', columns='Segment', values=['MeanDepth', 'percent_coverage']).reset_index()
    df_qc_pivot = df_qc_pivot.rename(columns=lambda x: f"Depth_{x}" if x != 'Sample' else x)
    # Pivot MeanDepth and percent_coverage columns
    df_qc_pivot = df_qc_cleaned.pivot(index='Sample', columns='Segment', values=['MeanDepth', 'percent_coverage']).reset_index()

    # Flatten hierarchical column names
    df_qc_pivot.columns = ['_'.join(col).strip('_') if isinstance(col, tuple) else col for col in df_qc_pivot.columns]

    # Rename columns to have consistent names
    df_qc_pivot = df_qc_pivot.rename(columns=lambda x: x.replace('MeanDepth_', 'Depth_').replace('percent_coverage_', 'percent_coverage_'))

    df_qc_pivot = collapse_depth_columns(df_qc_pivot)

    # Read TSV file
    df_tsv = pd.read_csv(tsv_file_path, sep='\t', usecols=['seqName', 'clade', 'qc.overallStatus', 'short-clade'])
    df_tsv_filtered = df_tsv[~df_tsv['clade'].isin(['NA', 'NaN', 'na', None])].dropna(subset=['clade'])
    df_tsv_filtered['Sample'] = df_tsv_filtered['seqName'].str.split('_').str[0]


    # Merge QC and TSV data
    df_merged = pd.merge(df_qc_pivot, df_tsv_filtered, on='Sample', how='left')

    # Read and merge Excel data
    df_excel = pd.read_excel(excel_file_path, sheet_name='1_Subtype Predictions')
    df_excel = df_excel.drop(columns=[col for col in df_excel.columns if "H: type prediction" in col or "N: type prediction" in col])
    df_excel['Type'] = df_excel['H: top match virus name'].apply(
        lambda x: 'Influenza A virus' if pd.notnull(x) and 'Influenza A' in x else (
            'Influenza B virus' if pd.notnull(x) and 'Influenza B' in x else 'Other'
        )
    )
    df_excel.rename(columns={'Subtype Prediction': 'Subtype'}, inplace=True)
    df_excel['Sample'] = df_excel['Sample'].astype(str)
    df_merged = pd.merge(df_merged, df_excel, on='Sample', how='left')

    # Read and merge Kraken2 data
    df_kraken2 = pd.read_csv(kraken2_file_path, sep='\t')
    df_kraken2['Sample'] = df_kraken2['Sample'].astype(str)
    df_merged = pd.merge(df_merged, df_kraken2, on='Sample', how='left')

    # Read and merge drug sensitivity data
    with open(drug_sensitivity_path, 'r') as file:
        lines = file.readlines()
    drug_data_end_index = next(i for i, line in enumerate(lines) if line.startswith('Counts of Individual Drug Sensitivity by subtype'))
    df_drug = pd.read_csv(io.StringIO('\n'.join(lines[:drug_data_end_index])), sep='\t')
    df_drug['Sample'] = df_drug['Query'].str.extract(r'(\d+)')[0]
    df_drug_pivot = df_drug.pivot_table(index='Sample', columns='drugname', values='sensitivity', aggfunc='first')
    df_drug_pivot.columns = [f"{col}_sensitivity" for col in df_drug_pivot.columns]
    df_merged = pd.merge(df_merged, df_drug_pivot.reset_index(), on='Sample', how='left')

    # Read and merge mutations
    df_mutations = aggregate_mutations(mut)
    print(df_mutations.head())
    df_merged = pd.merge(df_merged, df_mutations, on='Sample', how='left')

    # Read and merge flu typing report
    df_typing = pd.read_csv(typing_report_path, sep='\t')
    df_typing['Sample'] = df_typing['Sample'].astype(str)
    df_merged = pd.merge(df_merged, df_typing, on='Sample', how='left')

    # Read and merge IRMA QC data
    df_irma = pd.read_csv(irma_qc_path, sep='\t')
    df_irma['Sample'] = df_irma['Sample'].astype(str)
    df_merged = pd.merge(df_merged, df_irma, on='Sample', how='left')
    
    # Add RunID and RunDate
    df_merged['RunID'] = f'{run_id}_IAV'
    df_merged['RunDate'] = datetime.now().strftime('%Y-%m-%d')

    # Genomic QC flag
    df_merged['Genomic_QC'] = df_merged.apply(
        lambda row: 'PASS' if row.get('Depth_HA', 0) > 30 and row.get('Depth_NA', 0) > 30 else 'FAIL', axis=1
    )


    # Reorder columns
    desired_order = ['Sample', 'RunID', 'RunDate'] + \
    [col for col in df_merged.columns if 'Depth_' in col or 'percent_coverage_' in col] + \
    ['Genomic_QC'] + \
    [col for col in df_irma.columns if col != 'Sample'] + \
    [col for col in df_kraken2.columns if col != 'Sample'] + \
    [col for col in df_typing.columns if col != 'Sample'] + \
    [col for col in df_tsv_filtered.columns if col != 'Sample'] + \
    [col for col in df_excel.columns if col != 'Sample'] + \
    [col for col in df_drug_pivot.columns if col != 'Sample'] + \
    [col for col in df_mutations.columns if col != 'Sample']

    desired_order = list(dict.fromkeys(desired_order))
    df_merged = df_merged[[col for col in desired_order if col in df_merged.columns]]

    # Export to CSV
    output_csv_path = f'{run_id}.csv'
    df_merged.to_csv(output_csv_path, index=False)
    
    # Import to SQL
    import_to_sql(df_merged, 'resp_results')


    return output_csv_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Import TSV and Excel data into an MS SQL database table.')
    parser.add_argument('--tsv', required=True, help='Path to the TSV file.')
    parser.add_argument('--excel', required=True, help='Path to the Excel file.')
    parser.add_argument('--qc', required=True, help='Path to the QC file.')
    parser.add_argument('--kraken2', required=True, help='Path to the Kraken2 file.')
    parser.add_argument('--drug', required=True, help='Path to the drug sensitivity TSV file.')
    parser.add_argument('--mut', required=True, help='Path to the output text file.')
    parser.add_argument('--typing', required=True, help='Path to the typing report.')
    parser.add_argument('--run', required=True, help='RunID must be entered e.g. VON23003')
    parser.add_argument('--irma', required=True, help='Path to irma QC file.')

    args = parser.parse_args()
    output_csv_path = main(args.tsv, args.excel, args.qc, args.kraken2, args.drug, args.typing, args.irma, args.run, args.mut)

