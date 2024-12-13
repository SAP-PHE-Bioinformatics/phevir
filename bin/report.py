#!/usr/bin/env python3
import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from argparse import ArgumentParser

sns.set_style('darkgrid')

def create_report(df_sample, df_resistance, output_file):
    with PdfPages(output_file) as pdf:
        fig, ax = plt.subplots(figsize=(8.5, 11))
        ax.axis('off')

        # Title
        fig.suptitle('Influenza Workflow Report', fontsize=16)

        # Sample Details Table
        ax.text(0.5, 0.9, 'Sample Details', ha='center', fontsize=14)
        table_sample = ax.table(cellText=df_sample.values, colLabels=df_sample.columns, cellLoc='center', loc='center', bbox=[0.1, 0.6, 0.8, 0.3])
        table_sample.auto_set_font_size(False)
        table_sample.set_fontsize(10)
        
        # Resistance Details Table
        ax.text(0.5, 0.45, 'Resistance Details', ha='center', fontsize=14)
        table_resistance = ax.table(cellText=df_resistance.values, colLabels=df_resistance.columns, cellLoc='center', loc='center', bbox=[0.1, 0.05, 0.8, 0.3])
        table_resistance.auto_set_font_size(False)
        table_resistance.set_fontsize(10)

        pdf.savefig()
        plt.close()

def run(opts):
    # Load input data
    df_subtyping = pd.read_csv(opts.subtyping)
    df_mosdepth = pd.read_csv(opts.mosdepth)
    df_minimap2 = pd.read_csv(opts.minimap2_stats)
    df_bcftools = pd.read_csv(opts.bcftools_stats)
    df_nextclade = pd.read_csv(opts.nextclade_results)

    # Merge dataframes to create the sample details table
    df_sample = df_subtyping.merge(df_mosdepth, on='SampleID')
    df_sample = df_sample.merge(df_minimap2, on='SampleID')
    df_sample = df_sample.merge(df_bcftools, on='SampleID')
    df_sample = df_sample.merge(df_nextclade, on='SampleID')

    df_sample['DOB'] = pd.to_datetime(df_sample['DOB']).dt.strftime('%Y-%m-%d')

    sample_columns = ['SampleID', 'patientName', 'DOB', 'Influenza Type', 'Subtyping', 'nextclade_clade']
    df_sample = df_sample[sample_columns]

    # Create the resistance details table
    df_resistance = df_bcftools[['SampleID', 'subtype', 'Mutation', 'invivo/invitro', 'Oseltamivir', 'Zanamivir', 'Peramivir', 'Chromosome gene']]

    # Create PDF report
    create_report(df_sample, df_resistance, opts.output_pdf)

if __name__ == '__main__':
    parser = ArgumentParser(description='Create Influenza Workflow Report')
    parser.add_argument('-st', '--subtyping', required=True, help='Subtyping report CSV')
    parser.add_argument('-md', '--mosdepth', required=True, help='Mosdepth results CSV')
    parser.add_argument('-mm2', '--minimap2_stats', required=True, help='Minimap2 stats CSV')
    parser.add_argument('-bc', '--bcftools_stats', required=True, help='BCFTools stats CSV')
    parser.add_argument('-nc', '--nextclade_results', required=True, help='Nextclade results CSV')
    parser.add_argument('-o', '--output_pdf', required=True, help='Output PDF file')
    opts = parser.parse_args()

    run(opts)