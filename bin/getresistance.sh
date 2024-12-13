#!/usr/bin/env bash

# Define the file path
FASTA_FILE="$1"

# Verify the file exists
if [ ! -f "$FASTA_FILE" ]; then
  echo "FASTA file not found: $FASTA_FILE"
  exit 1
fi

# Submit the FASTA file and capture the response in a file
response_file="response.html"
curl -X POST \
  -H "Content-Type: multipart/form-data" \
  -F "seqfile=@${FASTA_FILE}" \
  -F "forceref=auto" \
  -F "lclq=1" \
  https://flusurver.bii.a-star.edu.sg/cgi-bin/flumapBlast3.pl -o $response_file

# Print the first few lines of the response file for debugging
echo "Response captured in $response_file"
head -n 20 $response_file

# Ensure the URLs are complete by appending the base URL
base_url="https://flusurver.bii.a-star.edu.sg"
mutation_report_url="${base_url}/${mutation_report_url#\.\./}"
query_summary_report_url="${base_url}/${query_summary_report_url#\.\./}"
query_clade_report_url="${base_url}/${query_clade_report_url#\.\./}"
drug_resistance_url="${base_url}/${drug_resistance_url#\.\./}"

# Step 2: Parse the response for the results URLs
mutation_report_url=$(grep -oP "(?<=href=')\.\./tmp/[^']+_result[^']+\.txt" $response_file)
query_summary_report_url=$(grep -oP "(?<=href=')\.\./tmp/[^']+_perquery\.csv" $response_file)
query_clade_report_url=$(grep -oP "(?<=href=')\.\./tmp/[^']+_perquery\.tsv" $response_file)
drug_sensitivity_report_url=$(grep -oP "(?<=href=')\.\./tmp/[^']+_drugsensitivity_summary[^']+\.tsv" $response_file)

# Ensure the URLs are complete by appending the base URL
base_url="https://flusurver.bii.a-star.edu.sg"
mutation_report_url="${base_url}/${mutation_report_url#\.\./}"
query_summary_report_url="${base_url}/${query_summary_report_url#\.\./}"
query_clade_report_url="${base_url}/${query_clade_report_url#\.\./}"
drug_sensitivity_report_url="${base_url}/${drug_sensitivity_report_url#\.\./}"

# Print the parsed URLs for debugging
echo "Mutation Report URL: $mutation_report_url"
echo "Query Summary Report URL: $query_summary_report_url"
echo "Query Clade Report URL: $query_clade_report_url"
echo "Drug Sensitivity Report URL: $drug_sensitivity_report_url"

# Step 3: Download the results if URLs are found
if [ -n "$mutation_report_url" ]; then
  curl -X GET "$mutation_report_url" -o mutation_report.txt
  echo "Mutation report downloaded: mutation_report.txt"
else
  echo "Mutation report URL not found in response."
fi

if [ -n "$query_summary_report_url" ]; then
  curl -X GET "$query_summary_report_url" -o query_summary_report.csv
  echo "Query summary report downloaded: query_summary_report.csv"
else
  echo "Query summary report URL not found in response."
fi

if [ -n "$query_clade_report_url" ]; then
  curl -X GET "$query_clade_report_url" -o query_clade_report.tsv
  echo "Query clade report downloaded: query_clade_report.tsv"
else
  echo "Query clade report URL not found in response."
fi

if [ -n "$drug_sensitivity_report_url" ]; then
  curl -X GET "$drug_sensitivity_report_url" -o drug_sensitivity_report.tsv
  echo "Drug sensitivity report downloaded: drug_sensitivity_report.tsv"
else
  echo "Drug sensitivity report URL not found in response."
fi

echo "Results downloaded."