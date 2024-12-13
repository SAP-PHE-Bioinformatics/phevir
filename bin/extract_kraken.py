import os
import sys
import pandas as pd

def get_top_matches(directory, output_tsv):
    # Check if the directory exists
    if not os.path.isdir(directory):
        print(f"The directory {directory} does not exist.")
        return
    
    # List to store results
    results = []
    
    # Walk through the directory to find .report files
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".report.txt"):
                file_path = os.path.join(root, file)
                sample_name = os.path.splitext(file)[0]
                top_matches = extract_top_matches(file_path)
                results.append([sample_name] + top_matches)
    
    # Create DataFrame and save to TSV
    headers = ["Sample", "1st Match", "2nd Match", "3rd Match"]
    df = pd.DataFrame(results, columns=headers)
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Results saved to {output_tsv}")

def extract_top_matches(file_path):
    # Read the report file into a DataFrame
    df = pd.read_csv(file_path, sep='\t', header=None, names=[
        'percentage', 'reads_clade', 'reads_direct', 'taxonomy_rank', 'taxonomy_id', 'taxonomy_name'
    ])
    
    # Filter for genus level
    genus_df = df[df['taxonomy_rank'] == 'G']
    
    # Sort by percentage and get the top 3 matches
    top_matches = genus_df.nlargest(3, 'percentage')[['taxonomy_name', 'percentage']]
    
    # Format the top matches as strings and bold 'Alphainfluenzavirus'
    top_matches_list = [
        f"**{row['taxonomy_name']} ({row['percentage']}%)**" if row['taxonomy_name'] == 'Alphainfluenzavirus' else f"{row['taxonomy_name']} ({row['percentage']}%)"
        for _, row in top_matches.iterrows()
    ]
    
    # Ensure the list has exactly 3 elements
    while len(top_matches_list) < 3:
        top_matches_list.append("")
    
    # Strip any leading/trailing whitespace from elements
    top_matches_list = [match.strip() for match in top_matches_list]
    
    return top_matches_list

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <directory> <output_tsv>")
        sys.exit(1)
    
    directory = sys.argv[1]
    output_tsv = sys.argv[2]
    get_top_matches(directory, output_tsv)
