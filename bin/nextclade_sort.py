#!/usr/bin/env python

import os
import csv
import argparse

# Dictionary containing data for various flu subtypes with their associated datasets.
flu_subtypes = {
    "H1N1": {"dataset": "flu_h1n1pdm_ha"},
    "H3N2": {"dataset": "flu_h3n2_ha"},
    "Victoria": {"dataset": "flu_vic_ha"},
    "Yamagata": {"dataset": "flu_yam_ha"},
    "H5N1": {"dataset": "community/moncla-lab/iav-h5/ha/all-clades"}
}


def main():
    # Set up an argument parser to accept input arguments for the script
    parser = argparse.ArgumentParser(
        description="Outputs the dataset, reference, and tag for the HA gene of a given flu subtype."
    )
    parser.add_argument("--sample", required=True, help="Sample name")
    args = parser.parse_args()

    # Construct the path of the input file based on the provided sample name
    input_file_path = f"{args.sample}.combined.typing.tsv"

    # Check if the input file exists
    if not os.path.exists(input_file_path):
        print(f"Error: Input file '{input_file_path}' does not exist")
        return

    # Read the input file and find the flu subtype
    flu_subtype = None
    with open(input_file_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["Sample"] == args.sample:
                # Check for valid subtype in 'abricate_InsaFlu_subtype' and fallback to 'IRMA_subtype' if needed
                abricate_subtype = row["abricate_InsaFlu_subtype"].strip()
                irma_subtype = row["IRMA_subtype"].strip()

                # Handle empty or "No abricate subtype" in the 'abricate_InsaFlu_subtype' column
                if abricate_subtype and abricate_subtype.lower() != "no abricate subtype":
                    flu_subtype = abricate_subtype
                else:
                    flu_subtype = irma_subtype
                break

    # If no flu subtype is found, handle it as an error
    if not flu_subtype:
        print(f"Error: No valid subtype found for sample '{args.sample}'")
        return

    # Verify that the found subtype is valid and present in the flu_subtypes dictionary
    if flu_subtype not in flu_subtypes:
        print(f"Error: Invalid flu subtype '{flu_subtype}' for sample '{args.sample}'")
        return

    # Prepare the dataset and output file path
    dataset = flu_subtypes[flu_subtype]["dataset"]
    output_file_path = f"{args.sample}_dataset.txt"

    # Write to the output file and print information
    with open(output_file_path, "w") as f:
        f.write(f"{dataset}\n")
        print(f"  {dataset}: {dataset} (output to {output_file_path})")


if __name__ == "__main__":
    main()
