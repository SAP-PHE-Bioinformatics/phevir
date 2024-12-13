#!/bin/python

import pandas as pd
import glob
import os
import sys
from os.path import exists

# File paths from user input
file_paths = [
    '/scratch/viral/nf-core-phevir/work/74/f975af8773d5ca2dbf56f6e3070387/2416914404.tsv',
    '/scratch/viral/nf-core-phevir/work/9a/928fbd32140da17bbca9644f6cb980/2416914404.stats',
    '/scratch/viral/nf-core-phevir/work/c0/1cc0dd6daaf23b9f6fadd8565de9af/2416512523.stats',
    '/scratch/viral/nf-core-phevir/work/18/ca77aa9153fd86d527250fff777416/2416512523.tsv',
    # (other file paths here)
]

# Constants
min_depth = int(sys.argv[1])

# Initialize dataframes
summary_df = pd.DataFrame(columns=['sample_id'])
columns = []

# Function to process FASTA files
def process_fasta_files(file_paths):
    fasta_df = pd.DataFrame(columns=['fasta_sample', "fasta_line", "num_N", "num_total"])
    for file in file_paths:
        if file.endswith(('.fa', '.fasta', '.fna')):
            print("Getting basic information from fasta " + file)
            if os.stat(file).st_size != 0:
                sample = str(file).replace(".fasta", '').replace(".fna", '').replace('.fa', '').replace('.consensus', '')
                with open(file) as fasta:
                    fasta_line = ''
                    sequence = ''
                    for line in fasta:
                        if ">" in line:
                            fasta_line = line.replace(">", "").strip().replace(".consensus_threshold*", "")
                        else:
                            sequence += line.strip()
                    num_N = sequence.count('N') + sequence.count('n')
                    num_total = len(sequence)
                    tmp_fasta_df = pd.DataFrame({'fasta_sample': [sample], 'fasta_line': [fasta_line], 'num_N': [num_N], 'num_total': [num_total]})
                    fasta_df = pd.concat([fasta_df, tmp_fasta_df], axis=0)
    return fasta_df

# Process and merge FASTA data
fasta_df = process_fasta_files(file_paths)
if not fasta_df.empty:
    summary_df = pd.merge(summary_df, fasta_df, left_on='sample_id', right_on='fasta_sample', how='outer')
    summary_df['sample_id'] = summary_df['sample_id'].fillna(summary_df['fasta_sample'])
    summary_df = summary_df.drop('fasta_sample', axis=1)
    columns = ['fasta_line'] + columns + ['num_N', 'num_total']

summary_df['sample_id'] = summary_df['sample_id'].astype(object)

# Function to process KRAKEN2 files
def process_kraken2_files(file_paths):
    kraken2_df = pd.DataFrame(columns=['kraken2_sample', '%_human_reads', 'top_organism', 'percent_reads_top_organism'])
    for file in file_paths:
        if file.endswith('_kraken2_report.txt'):
            print("Getting species information from " + file)
            if os.stat(file).st_size != 0:
                sample = str(file).replace('_kraken2_report.txt', '')
                percent_human_reads = 0
                top_organism = 'none'
                percent_reads_top_organism = 0
                with open(file) as report:
                    kraken2_sample_df = pd.DataFrame(columns=['percent', 'reads', 'species'])
                    for line in report:
                        if "S" == line.split()[3]:
                            per = line.split()[0].strip()
                            reads = int(line.split()[1].strip())
                            species = ' '.join(line.split()[5:]).strip()
                            tmp_kraken2_df = pd.DataFrame({'percent': [per], 'reads': [reads], 'species': [species]})
                            kraken2_sample_df = pd.concat([kraken2_sample_df, tmp_kraken2_df], axis=0)
                    kraken2_sample_df = kraken2_sample_df.sort_values(by=['reads'], ascending=False)
                    kraken2_human_df = kraken2_sample_df[kraken2_sample_df['species'].isin(['Homo sapiens', 'Human'])]
                    if not kraken2_human_df.empty:
                        percent_human_reads = kraken2_human_df['percent'].iloc[0]
                    kraken2_species_df = kraken2_sample_df[~kraken2_sample_df['species'].isin(['Homo sapiens', 'Human'])]
                    if not kraken2_species_df.empty:
                        top_organism = kraken2_species_df['species'].iloc[0]
                        percent_reads_top_organism = kraken2_species_df['percent'].iloc[0]
                    kraken2_tmp_df = pd.DataFrame({'kraken2_sample': [sample], '%_human_reads': [percent_human_reads], 'top_organism': [top_organism], 'percent_reads_top_organism': [percent_reads_top_organism]})
                    kraken2_df = pd.concat([kraken2_df, kraken2_tmp_df], axis=0)
    return kraken2_df

# Process and merge KRAKEN2 data
kraken2_df = process_kraken2_files(file_paths)
if not kraken2_df.empty:
    summary_df = pd.merge(summary_df, kraken2_df, left_on='sample_id', right_on='kraken2_sample', how='outer')
    summary_df['sample_id'] = summary_df['sample_id'].fillna(summary_df['kraken2_sample'])
    summary_df = summary_df.drop('kraken2_sample', axis=1)
    columns = columns + ['top_organism', 'percent_reads_top_organism', '%_human_reads']

# Function to process DEPTH files
def process_depth_files(file_paths, min_depth):
    depth_df = pd.DataFrame(columns=['samtools_sample', "num_pos_" + str(min_depth) + "X"])
    for file in file_paths:
        if file.endswith('.depth.txt'):
            print("Finding bases above " + str(min_depth) + " in " + file)
            if os.stat(file).st_size != 0:
                ind_depth_df = pd.read_table(file, header=None)
                depth = ind_depth_df[2][ind_depth_df[2] > min_depth].count()
                sample = str(file).replace('.depth.txt', '')
                tmp_depth_df = pd.DataFrame({'samtools_sample': [sample], "num_pos_" + str(min_depth) + "X": [depth]})
                depth_df = pd.concat([depth_df, tmp_depth_df], axis=0)
    return depth_df

# Process and merge DEPTH data
depth_df = process_depth_files(file_paths, min_depth)
if not depth_df.empty:
    summary_df = pd.merge(summary_df, depth_df, left_on='sample_id', right_on='samtools_sample', how='outer')
    summary_df['sample_id'] = summary_df['sample_id'].fillna(summary_df['samtools_sample'])
    summary_df = summary_df.drop('samtools_sample', axis=1)
    depth_columns = list(depth_df.columns)
    depth_columns.remove('samtools_sample')
    columns = columns + depth_columns

# Process ACI file
aci_file = next((f for f in file_paths if f.endswith('aci_coverage_summary.csv')), None)
if aci_file and exists(aci_file):
    print("Getting results from ACI " + aci_file)
    aci_df = pd.read_csv(aci_file, dtype=str)
    aci_df = aci_df.add_prefix('aci_')
    tmp_df = aci_df.iloc[:, aci_df.columns != 'aci_bam'].astype('float').copy()
    aci_df['aci_num_failed_amplicons'] = tmp_df.apply(lambda x: x[x < min_depth].count(), axis=1)
    aci_df['aci_name'] = aci_df['aci_bam'].str.replace(".primertrim.sorted.bam", "", regex=False).str.replace('.sorted.bam', '', regex=False)
    summary_df = pd.merge(summary_df, aci_df[['aci_name', 'aci_num_failed_amplicons']], left_on='sample_id', right_on='aci_name', how='outer')
    summary_df['sample_id'].fillna(summary_df['aci_name'], inplace=True)
    summary_df.drop('aci_name', axis=1, inplace=True)
    columns = columns + ['aci_num_failed_amplicons']

# Process Samtools Ampliconstats file
ampliconstats_file = next((f for f in file_paths if f.endswith('ampliconstats.summary')), None)
if ampliconstats_file and exists(ampliconstats_file):
    print("Getting results from samtools ampliconstats file " + ampliconstats_file)
    amp_df = pd.read_table(ampliconstats_file, header=None)
    amp_df = amp_df.rename(columns={1: 'sample'})
    amp_df['samtools_num_failed_amplicons'] = amp_df.iloc[:, 2:].lt(20).sum(axis=1)
    amp_df['sample_name'] = amp_df['sample'].str.replace('.primertrim.sorted', '', regex=False)
    amp_tmp_df = amp_df[['sample_name', 'samtools_num_failed_amplicons']]
    summary_df = pd.merge(summary_df, amp_tmp_df, left_on='sample_id', right_on='sample_name', how='outer')
    summary_df['sample_id'].fillna(summary_df['sample_name'], inplace=True)
    summary_df.drop('sample_name', axis=1, inplace=True)
    columns = columns + ['samtools_num_failed_amplicons']

# Process QUAST files
def process_quast_files(file_paths):
    quast_df = pd.DataFrame(columns=['quast_sample', 'N50', 'num_contigs'])
    for file in file_paths:
        if file.endswith('.tsv'):
            print("Getting QUAST results from " + file)
            if os.stat(file).st_size != 0:
                df = pd.read_csv(file, sep='\t')
                # Assuming the file has columns 'N50', 'num_contigs'
                sample = os.path.basename(file).replace('.tsv', '')
                N50 = df['N50'].iloc[0] if 'N50' in df.columns else None
                num_contigs = df['num_contigs'].iloc[0] if 'num_contigs' in df.columns else None
                tmp_quast_df = pd.DataFrame({'quast_sample': [sample], 'N50': [N50], 'num_contigs': [num_contigs]})
                quast_df = pd.concat([quast_df, tmp_quast_df], axis=0)
    return quast_df

ivar_variants_df = pd.DataFrame(columns=['ivar_sample', 'ivar_num_variants_identified'])
ivar_variant_files = glob.glob("*.variants.tsv")
for file in ivar_variant_files :
    print("Counting variants in " + file)
    if not ( os.stat(file).st_size==0 ) :
        ivar_variants = 0
        lines = open(file, 'r').readlines()
        for line in lines :
            if "TRUE" in line.split("\t")[13]:
                ivar_variants += 1

        sample              = str(file).replace('.variants.tsv', '')
        tmp_ivar_df         = pd.DataFrame({'ivar_sample': [sample], 'ivar_num_variants_identified': [ivar_variants]})
        ivar_variants_df    = pd.concat([ivar_variants_df, tmp_ivar_df], axis=0 )

if not ivar_variants_df.empty :
    summary_df              = pd.merge(summary_df, ivar_variants_df, left_on = 'sample_id', right_on = 'ivar_sample', how = 'outer')
    summary_df['sample_id'] = summary_df['sample_id'].fillna(summary_df['ivar_sample'])
    summary_df              = summary_df.drop('ivar_sample', axis=1)
    ivar_variants_columns   = ['ivar_num_variants_identified']
    columns                 = columns + ivar_variants_columns


# Process and merge QUAST data
quast_df = process_quast_files(file_paths)
if not quast_df.empty:
    summary_df = pd.merge(summary_df, quast_df, left_on='sample_id', right_on='quast_sample', how='outer')
    summary_df['sample_id'] = summary_df['sample_id'].fillna(summary_df['quast_sample'])
    summary_df = summary_df.drop('quast_sample', axis=1)
    columns = columns + ['N50']

def vadr_sample_name(s):
    if s.count('.') >=1:
        if len(s.split(".")[-1]) > 2:
            return ''.join(s.split(".")[:-1])
    return s

if exists(vadr_file) :
    print("Getting results from vadr file " + vadr_file)
    vadr_df = pd.read_csv(vadr_file, dtype = str, usecols = ['name', 'p/f', 'model', 'alerts'], index_col= False)
    vadr_df = vadr_df[vadr_df['name'] != 'name']
    vadr_df = vadr_df[vadr_df['name'] != 'seq']
    vadr_df = vadr_df.add_prefix('vadr_')
    vadr_columns = list(vadr_df.columns)
    vadr_columns.remove('vadr_name')
    vadr_columns.remove('vadr_p/f')

    if 'fasta_line' in summary_df.columns.tolist():
        summary_df = pd.merge(summary_df, vadr_df, left_on = 'fasta_line', right_on = 'vadr_name', how = 'outer')
        summary_df['sample_id'].fillna(summary_df['vadr_name'], inplace=True)
        summary_df.drop('vadr_name', axis=1, inplace=True)
        columns = ['vadr_p/f'] + columns + vadr_columns
    else:
        vadr_df['sample_match'] = vadr_df['vadr_name'].str.replace('Consensus_', '', regex =  False).apply(vadr_sample_name)

        summary_df = pd.merge(summary_df, vadr_df, left_on = 'sample_id', right_on = 'sample_match', how = 'outer')
        summary_df['sample_id'].fillna(summary_df['sample_match'], inplace=True)
        summary_df.drop('vadr_name', axis=1, inplace=True)
        summary_df.drop('sample_match', axis=1, inplace=True)
        columns = ['vadr_p/f'] + columns + vadr_columns 

if exists(nextclade_file) :
    print("Getting results from nextclade file " + nextclade_file)

    use_cols = ['seqName', 'clade', 'qc.overallStatus', 'qc.overallScore']

    first = pd.read_table(nextclade_file, sep = ';' , dtype = str, nrows=1)
    if 'clade_who' in first.columns:
        use_cols.append('clade_who')
    if 'outbreak' in first.columns:
        use_cols.append('outbreak')
    if 'lineage' in first.columns:
        use_cols.append('lineage')

    nextclade_df = pd.read_table(nextclade_file, sep = ';' , dtype = str, usecols = use_cols)
    nextclade_df=nextclade_df.add_prefix('nextclade_')
    nextclade_columns = list(nextclade_df.columns)
    nextclade_df['sample_match'] = nextclade_df['nextclade_seqName'].str.replace('Consensus_', '', regex =  False).str.split(' ').str[0]
    nextclade_columns.remove('nextclade_seqName')
    nextclade_columns.remove('nextclade_clade')

    summary_df = pd.merge(summary_df, nextclade_df, left_on = 'sample_id', right_on = 'sample_match', how = 'outer')
    summary_df['sample_id'].fillna(summary_df['sample_match'], inplace=True)
    summary_df.drop('nextclade_seqName', axis=1, inplace=True)
    summary_df.drop('sample_match', axis = 1, inplace = True )
    columns = ['nextclade_clade'] + columns + nextclade_columns


if exists(versions_file) :
    print("Adding versions to summary file from " + versions_file)
    versions_df = pd.read_csv(versions_file, dtype = str)
    version_columns = list(versions_df.columns)

    for version in version_columns:
        software_version = versions_df[version].iloc[0]
        summary_df[version] = software_version
    columns = columns + version_columns

summary_df['sample'] = summary_df['sample_id'].str.split("_").str[0]
summary_df['sample_id'] = summary_df['sample_id'].astype('string')
summary_df = summary_df.replace([" ", ",", "\t", "\n"], [" ", " ", " ", " "], regex=True)
summary_df = summary_df.sort_values(by=['sample_id'], ascending=True)
summary_df = summary_df.replace([" ", ",", "\t", "\n"], [" ", " ", " ", " "], regex=True)
summary_df = summary_df.drop_duplicates(keep='first')
summary_df.to_csv('results.csv', columns=['sample_id', 'sample'] + columns, index=False)
summary_df.to_csv('results.txt', columns=['sample_id', 'sample'] + columns, index=False, sep="\t")
