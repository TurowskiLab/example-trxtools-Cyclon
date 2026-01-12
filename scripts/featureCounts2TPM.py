import pandas as pd
import os
import sys
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Convert featureCounts output to TPM.')
parser.add_argument('-i', '--input', required=True, help='Path to the input featureCounts file.')
parser.add_argument('-o', '--output', required=True, help='Path to the output TPM file.')
args = parser.parse_args()

input_file = args.input
output_file = args.output

if not os.path.isfile(input_file):
    print(f"File not found: {input_file}")
    sys.exit(1)
else:
    # Read and wrangle featureCounts output
    df = pd.read_csv(input_file, sep="\t", skiprows=1)
    df = df.drop(columns=['Chr', 'Start', 'End', 'Strand']).set_index('Geneid')
    # Calculate TPM
    df.iloc[:, 1:] = df.iloc[:, 1:].div(df['Length']/1e3, axis=0)
    df.iloc[:, 1:] = df.iloc[:, 1:].div(df.iloc[:, 1:].sum()/1e6, axis=1)
    df = df.drop(columns=['Length'])
    df.to_csv(output_file, sep="\t")
