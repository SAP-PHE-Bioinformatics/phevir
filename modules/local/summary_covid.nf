process SUMMARY {
    tag "${batch_id}-generate_summary"
    label 'process_medium'

    input:
    path(QC) // Input coverage stats
    path(pangolin) // Pangolin results
    path(nextclade) // Nextclade results
    path(pangoCollapse) // Pango-collapse output

    output:
    path("${batch_id}_summary.tsv") // Raw TSV output
    path("${batch_id}_summary.xlsx") // Formatted Excel output

    script:
    """
    python -c "
import pandas as pd
import os

# Load input files
fastq_df = pd.read_csv('${QC}', sep='\\t', index_col=0)
pangolin_df = pd.read_csv('${pangolin}', sep='\\t', index_col=0)
nextclade_df = pd.read_csv('${nextclade}', sep='\\t', index_col=0)
pango_collapse_df = pd.read_csv('${pangoCollapse}', sep='\\t', index_col=0)

# Merge DataFrames
summary_df = bam_df.merge(fastq_df, left_index=True, right_index=True, how='outer')
summary_df = summary_df.merge(pangolin_df, left_index=True, right_index=True, how='outer')
summary_df = summary_df.merge(nextclade_df, left_index=True, right_index=True, how='outer')
summary_df = summary_df.merge(pango_collapse_df, left_index=True, right_index=True, how='outer')

# Select and rename columns of interest
summary_columns = ['mean', 'median', 'Genome fraction (%)']
summary_df = summary_df[summary_columns].rename(columns={
    'mean': 'MeanDepthCov',
    'median': 'MedianDepthCov',
    'Genome fraction (%)': '%GenomeFrac'
}).fillna(0).sort_index()

summary_df['run_id'] = ${batch_id}
# Save as raw TSV
summary_df.to_csv('${batch_id}_COV_summary.tsv', sep='\\t', index=True, float_format='%.2f')

# Create formatted Excel output
with pd.ExcelWriter('${batch_id}_COV_summary.xlsx', engine='xlsxwriter') as writer:
    summary_df.to_excel(writer, index=True, header=True, sheet_name='Summary')
    workbook = writer.book
    worksheet = writer.sheets['Summary']

    # Define formats
    bad = workbook.add_format({'bold': False, 'font_color': '#E53935'}) # Red
    good = workbook.add_format({'bold': False, 'font_color': '#388E3C'}) # Green
    meh = workbook.add_format({'bold': False, 'font_color': '#ff9900'}) # Orange
    center = workbook.add_format({'align': 'center'})
    left = workbook.add_format({'align': 'left'})
    vertical_center = workbook.add_format({'align': 'center', 'valign': 'vcenter'})

    # Adjust column widths and formatting
    worksheet.set_column('A:A', 20, left)
    worksheet.set_column('B:C', 20, vertical_center)

    # Conditional formatting
    num_rows = len(summary_df.index) + 1
    worksheet.conditional_format(f'B2:B{num_rows}', {'type': 'cell', 'criteria': '>=', 'value': 90, 'format': good})
    worksheet.conditional_format(f'B2:B{num_rows}', {'type': 'cell', 'criteria': 'between', 'minimum': 65, 'maximum': 90, 'format': meh})
    worksheet.conditional_format(f'B2:B{num_rows}', {'type': 'cell', 'criteria': '<', 'value': 65, 'format': bad})

    writer.save()
    "
    """
}
