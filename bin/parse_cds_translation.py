#!/usr/bin/env python
import argparse
import pandas as pd
import csv
import re



def parse_resistance_mutations(resistance_mutation_list_df):
    resistance_mutations = []
    for index, row in resistance_mutation_list_df.iterrows():
        mutation = row['Mutation']
        mutations = mutation.split('+')
        resistance_mutations.append(mutations)
    return resistance_mutations

def parse_translation_file(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        current_sequence = ""
        current_id = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = current_sequence
                current_id = line[1:]
                current_sequence = ""
            else:
                current_sequence += line
        if current_id:
            sequences[current_id] = current_sequence
    return sequences

def detect_resistance_single_mutation(mutation, sequence):
    detected = False
    note = ""
    position = int(mutation[1:-1])  
    expected_amino_acid = mutation[-1]
    
    if sequence[position -1 ] == expected_amino_acid:  # account for 0 based indexing of position
        detected = True
        note += f"Detected mutation: {mutation}\n"
    else:
        note += f"Mutation {mutation} not detected. Found {sequence[position-1]} at position {position}.\n"

    return mutation, detected, note.strip()

def detect_resistance_combination_mutation(mutations, sequence):
    detected = True 
    note = ""

    for mutation in mutations:
        single_mutation_result = detect_resistance_single_mutation(mutation, sequence)
        mutation_detected = single_mutation_result[1]
        mutation_note = single_mutation_result[2]

        if not mutation_detected:
            detected = False

        note += mutation_note + "\n"

    if detected:
        detected_mutations = [mutation for mutation in mutations if detect_resistance_single_mutation(mutation, sequence)[1]]
        note += f"All mutations detected: {' '.join(detected_mutations)}\n"
    else:
        undetected_mutations = [mutation for mutation in mutations if not detect_resistance_single_mutation(mutation, sequence)[1]]
        note += f"Mutations not detected: {' '.join(undetected_mutations)}\n"

    note = note.strip().replace('\n', ' ').replace('..', '')

    return mutations, detected, note.strip()

def main():
    parser = argparse.ArgumentParser(description='Parse Nextclade results based on list of resistance mutations and reference amino acid')
    parser.add_argument('--cds_translation', type=str, required=True, help='Path to nextclade cds translation file for gene of interest.')
    parser.add_argument('--resistance_mutation_list', type=str, required=True, help='Path to list of resistance mutations')
    parser.add_argument('--output', type=str, required=True, help='Path to output CSV file')
    parser.add_argument('--gene', type=str, required=True, help='Gene of interest')
    args = parser.parse_args()

    # Import data
    sequences = parse_translation_file(args.cds_translation)
    resistance_mutation_list_df = pd.read_csv(args.resistance_mutation_list)

    # Parse the list of resistance mutations from input
    resistance_mutations = parse_resistance_mutations(resistance_mutation_list_df)



    results_list = []

    # For each sequence in the translation file
    for seq_name, sequence in sequences.items():
        # Clean sequence name if necessary
        seq_name = seq_name.replace(' Human respiratory syncytial virus A isolate hRSV/A/England/397/2017, complete genome', '')

        # For each mutation in the list of resistance mutations
        for mutations in resistance_mutations:
            result = None
            if len(mutations) > 1:
                result = detect_resistance_combination_mutation(mutations, sequence)
            else:
                result = detect_resistance_single_mutation(mutations[0], sequence)
            
            mutation_name = "+".join(mutations)
            detected = result[1]
            note = result[2]
            gene = str(args.gene)

            results_list.append([seq_name, gene, mutation_name, 'Present' if detected else 'Absent', note])

    results_df = pd.DataFrame(results_list, columns=['seqName', 'Gene', 'Mutation', 'Detected', 'Note'])
    results_df.to_csv(args.output, index=False, quoting=csv.QUOTE_NONE, escapechar='\\')

if __name__ == "__main__":
    main()
