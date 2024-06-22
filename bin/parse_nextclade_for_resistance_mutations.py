#!/usr/bin/env python
import argparse
import pandas as pd
import csv

def parse_resistance_mutations(resistance_mutation_list_df):
    # Getting mutations in list for combination of mutations (list of lists)
    resistance_mutations = []
    for index, row in resistance_mutation_list_df.iterrows():
        mutation = row['Mutation']
        mutations = mutation.split('+')
        resistance_mutations.append(mutations)
    return resistance_mutations

def parse_nextclade(nextclade_df, gene):
    # Parse nextclade results, select mutations in gene of interest and remove gene from mutation info
    aa_substitutions_per_sample = {}
    for index, row in nextclade_df.iterrows():
        seq_name = row['seqName'].replace(' Human respiratory syncytial virus A isolate hRSV/A/England/397/2017, complete genome', '')
        aa_substitutions = row['aaSubstitutions'].split(',')
        filtered_aa_substitutions = [sub.split(':', 1)[1] for sub in aa_substitutions if sub.startswith(gene + ':')]
        aa_substitutions_per_sample[seq_name] = filtered_aa_substitutions
    return aa_substitutions_per_sample

def detect_resistance_single_mutation(mutation, aa_substitutions, ref_positions_df):
    detected = False
    note = ""

    # Get info from input about reference position at mutation site 
    ref_amino_acid_at_position = ref_positions_df[ref_positions_df['Mutation'] == mutation]
    
    # Check that mutation is in input
    if not ref_amino_acid_at_position.empty:
        # Get reference info information (what base at mutation position)
        result = ref_amino_acid_at_position.iloc[0]['Result']

        # If WT aa is the same, then we can directly check for mutation in nextclade output
        if result == 'Match WT' and mutation in aa_substitutions:
            detected = True
            note += f"Detected mutation: {mutation}\n"
        
        # Otherwise, if the reference already contains the mutation, check for other mutations at the position of the mutation
        elif result == 'Match Mutant':
            position = int(mutation[1:-1])

            # If there are any mutations at that position, then the sample does not contain this mutation             
            if any(sub[1:-1] == mutation[1:-1] for sub in aa_substitutions):
                detected = False
                detected_mutations = [sub for sub in aa_substitutions if sub[1:-1] == mutation[1:-1]]
                note += f"Reference contains {mutation} but mutation has been detected at position {position}: Detected Mutation: {', '.join(detected_mutations)}.\n"
            # Otherwise, if no mutation in the sample, then the sample contains this mutation as well    
            else:
                detected = True
                note += f"Reference contains {mutation} and no mutations have been detected at position {position} so sample contains this mutation.\n"

        # In the case where the reference amino acid is not the same as the WT or mutant AA at this position, check if the mutant mutation is present
        # Example: For the mutation N63S, check if a 63S mutation is present. If mutant amino acid is present at that position, but we have to
        # ignore the amino acid at the first position since the reference differs from the wild type
        elif result == 'Match Neither':
            position = int(mutation[1:-1])
            substitution = mutation[-1]
            exact_mutation = f'{position}{substitution}'
            
            # Check if the position and amino acid substitution is present
            if any(sub[1:] == mutation[1:] for sub in aa_substitutions):
                detected = True
                note += f"Reference contains a different amino acid than the wild type at position {position}. Resistance mutation is {mutation} and a {exact_mutation} mutation has occurred.\n"
            # In this case, the sample contains a different amino acid than the wild type or mutant amino acid
            else:
                detected = False
                ref_amino_acid = ref_amino_acid_at_position.iloc[0]['Reference_AA']
                note += f"This sample contains amino acid {ref_amino_acid} at position {position}. This amino acid is the same as the reference which contains a different amino acid than the wild type but does not contain the mutant amino acid that confers resistance.\n"

    # If mutation is not in aa_substitutions and not found in ref_positions_df
    elif ref_amino_acid_at_position.empty:
        note += f"ERROR: Mutation missing from --ref_amino_acid input. Please supply reference amino acid information to get mutation detection information.\n"

    return mutation, detected, note.strip()

def detect_resistance_combination_mutation(mutations, aa_substitutions, ref_positions_df):

    detected = True 
    note = ""

    # Determine result for each individual mutation
    for mutation in mutations:
        single_mutation_result = detect_resistance_single_mutation(mutation, aa_substitutions, ref_positions_df)
        mutation_detected = single_mutation_result[1]
        mutation_note = single_mutation_result[2]

        # If any mutation is not detected, then combo mutation is false
        if not mutation_detected:
            detected = False

        note += mutation_note + "\n"

    # Output which mutations detected out of combo
    if detected:
        detected_mutations = [mutation for mutation in mutations if detect_resistance_single_mutation(mutation, aa_substitutions, ref_positions_df)[1]]
        note += f"All mutations detected: {' '.join(detected_mutations)}\n"
    else:
        undetected_mutations = [mutation for mutation in mutations if not detect_resistance_single_mutation(mutation, aa_substitutions, ref_positions_df)[1]]
        note += f"Mutations not detected: {' '.join(undetected_mutations)}\n"

    note = note.strip().replace('\n', ' ').replace('..', '')

    return mutations, detected, note.strip()

def main():
    parser = argparse.ArgumentParser(description='Parse Nextclade results based on list of resistance mutations and reference amino acid')
    parser.add_argument('--ref_amino_acid', type=str, required=True, help='Path to the file that contains the information on the reference amino acid at mutation positions. Required Result column with whether AA "Match WT", "Match Mutant" or "Match Neither" ')
    parser.add_argument('--nextclade', type=str, required=True, help='Path to nextclade_qc.tsv containing aaSubstitutions')
    parser.add_argument('--resistance_mutation_list', type=str, required=True, help='Path to the list of resistance mutations')
    parser.add_argument('--output', type=str, required=True, help='Path to the final output CSV file')
    args = parser.parse_args()

    gene_of_interest = "F"  # Can modify into script argument if needed. This will be the gene relating to the resistance mutations of interest.

    # Import data
    ref_mutation_positions_df = pd.read_csv(args.ref_amino_acid)
    nextclade_df = pd.read_csv(args.nextclade, sep='\t')
    resistance_mutation_list_df = pd.read_csv(args.resistance_mutation_list)



    # Parse the list of resistance mutations from input
    resistance_mutations = parse_resistance_mutations(resistance_mutation_list_df)
    # Parse the amino acid mutations detected in data
    aa_substitutions_per_sample = parse_nextclade(nextclade_df, gene_of_interest)

    print(aa_substitutions_per_sample)

    results_list = []

    # For each sample in nextclade results
    for seq_name, aa_subs in aa_substitutions_per_sample.items():

        
        # For each mutation in the list of resistance mutations
        for mutations in resistance_mutations:
            result = None
            # Check if combination of mutations or single mutation
            if len(mutations) > 1:
                result = detect_resistance_combination_mutation(mutations, aa_subs, ref_mutation_positions_df)
            else:
                result = detect_resistance_single_mutation(mutations[0], aa_subs, ref_mutation_positions_df)
            
            # For output file, recombine combination of mutations into original format
            mutation_name = "+".join(mutations)  
            detected = result[1]
            note = result[2]
            gene = gene_of_interest

            results_list.append([seq_name, gene, mutation_name, 'Present' if detected else 'Absent', note])

    # Output df with the resistance mutation, whether it was detected in results and note about results.
    results_df = pd.DataFrame(results_list, columns=['seqName', 'Gene', 'Mutation', 'Detected', 'Note'])
    results_df.to_csv(args.output, index=False, quoting=csv.QUOTE_NONE, escapechar='\\')

if __name__ == "__main__":
    main()