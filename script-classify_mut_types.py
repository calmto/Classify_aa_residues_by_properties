"""
____________________________
Title:
Analysis of a protein's drug binding site's tolerance/sensitivity by deep mutational scanning (DMS) and classification of amino acid properties into 'mutation types'

____________________________
Description:
This script requires a hardcoded wild-type (WT) sequence and a CSV file in the current directory.
The CSV file should contain the DMS results combined with the drug resistance screening of all mutants with the format:
 ,compound,seq_type,Nham_aa,aa_seq,s,cscore,refined_class,sensres
1,anidulafungin,single,0.0,FLVLSLRDP,0.06276288377908285,1.0,WT-like,sensitive
2,anidulafungin,single,1.0,*LVLSLRDP,-0.2306350142092216,1.0,WT-like,sensitive
3,anidulafungin,single,1.0,ALVLSLRDP,1.7562999789995668,1.0,intermediary,resistant
4,caspofungin,single,0.0,FLVLSLRDP,0.08276288377908285,1.0,WT-like,sensitive

This script was developed as part of a research project on preventive drug design of drug-resistant pathogens. 
The analysis involves: 1) creating mutant strains of every residue in the protein's drug binding site through DMS, 
2) measure how the mutant strains react to the drug(s) compared to a WT strain (degree of resistance/tolerance to a drug), 
3) Classify resistant mutations based on the amino acid properties of the WT and mutant residue, 
4) aggregate and plot the dataframe to a heatmap for easy visualization of drug-specific resistance mechanism.

Useful  for analyzing large datasets of single point mutations in a WT sequence based on a measured feature (resistance/tolerance in this case) and aa properties.

For more details, refer to the corresponding scientific article: [Link to the article]

____________________________
Written by: Alexandre Torbey, PhD Candidate, Department of Health and Biotechnology, Armand Frappier Institute, Quebec (Canada).
Main supervisors: Prof. David Chatenet (INRS-IAF) and Prof. Patrick Lague (Laval University)
Co-authors:  Dr. Romain Durand and Prof. Christian Landry (Laval University)
Date: 11/june/2024

____________________________
Dependencies:
- numpy
- pandas
- matplotlib

To install the necessary dependencies, run the following command:
pip install numpy pandas matplotlib

Ran on Python 3.10.10

____________________________
License:
This project is licensed under the MIT License - see the LICENSE file for details.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8

# Wild-type sequence
wt_sequence = "FLVLSLRDP"

# Load the CSV file into a DataFrame
file_path = 'FKS1-HS1_singles.csv'
df = pd.read_csv(file_path)

# Adjusted amino acid properties
aa_properties = {
    'Aromatic': "FYW",
    'Aliphatic': "ILVGAT",
    'Positive': "RK",
    'Negative': "DE",
    'Neutral': "NQSC",
    'Big': "ILFTWYKRQE",
    'Small': "VTDN",
    'Tiny': "AGCS",
    'Proline': "P",
    'Methionine': "M",
    'Histidine': "H"
}

# Define specific mutation transitions of interest and their inverses
mutation_transitions = {
    ("Aliphatic", "Aromatic"): "Aliphatic-Aromatic",
    ("Aliphatic", "Positive"): "Aliphatic-Positive",
    ("Aliphatic", "Negative"): "Aliphatic-Negative",
    ("Aliphatic", "Neutral"): "Aliphatic-Neutral",
    ("Aromatic", "Aliphatic"): "Aromatic-Aliphatic",
    ("Aromatic", "Positive"): "Aromatic-Positive",
    ("Aromatic", "Negative"): "Aromatic-Negative",
    ("Aromatic", "Neutral"): "Aromatic-Neutral",
    ("Positive", "Aliphatic"): "Positive-Aliphatic",
    ("Positive", "Aromatic"): "Positive-Aromatic",
    ("Positive", "Negative"): "Positive-Negative",
    ("Positive", "Neutral"): "Positive-Neutral",
    ("Negative", "Aliphatic"): "Negative-Aliphatic",
    ("Negative", "Aromatic"): "Negative-Aromatic",
    ("Negative", "Positive"): "Negative-Positive",
    ("Negative", "Neutral"): "Negative-Neutral",
    ("Neutral", "Aliphatic"): "Neutral-Aliphatic",
    ("Neutral", "Aromatic"): "Neutral-Aromatic",
    ("Neutral", "Positive"): "Neutral-Positive",
    ("Neutral", "Negative"): "Neutral-Negative",
    ("Big", "Small"): "Big-Small",
    ("Small", "Big"): "Small-Big",
    ("Big", "Tiny"): "Big-Tiny",
    ("Tiny", "Big"): "Tiny-Big",
    ("Small", "Tiny"): "Small-Tiny",
    ("Tiny", "Small"): "Tiny-Small",
    ("Any", "Proline"): "Any-Proline",
    ("Any", "Methionine"): "Any-Methionine",
    ("Any", "Histidine"): "Any-Histidine",
    ("Proline", "Any"): "Proline-Any",
    ("Methionine", "Any"): "Methionine-Any",
    ("Histidine", "Any"): "Histidine-Any",
}
ordered_mutation_types = list(mutation_transitions.values())
ordered_mutation_types.append("Similar to WT")

# Filtering the DataFrame to include only resistant mutations
df_resistant = df[(df['seq_type'] == 'single') & 
                  (df['sensres'] == 'resistant') & 
                  (df['compound'] != 'none') & 
                  (~df['aa_seq'].str.contains('\*'))].copy()

# Function to find the mutation position (1-based indexing)
def find_mutation_position(mutant_seq):
    for i, (wt_aa, mut_aa) in enumerate(zip(wt_sequence, mutant_seq)):
        if wt_aa != mut_aa:
            return i + 1
    return None

# Function to convert compound names to initials
def compound_to_initials(compound):
    initials = {'anidulafungin': 'A', 'caspofungin': 'C', 'micafungin': 'M'}
    return ''.join(sorted({initials.get(c, '') for c in compound.split(',')}))

# Function that handles mutation type based on position to determine resistance
def determine_mutation_type_and_resistance(wt_seq, mutant_seq, position, compound):
    if position is None:
        return [('WT', compound_to_initials(compound))]
    wt_aa = wt_seq[position - 1]
    mut_aa = mutant_seq[position - 1]
    applicable_types = []
    wt_groups = set(group for group, aas in aa_properties.items() if wt_aa in aas)
    mut_groups = set(group for group, aas in aa_properties.items() if mut_aa in aas)
    primary_props = {'Aliphatic', 'Aromatic', 'Positive', 'Negative', 'Neutral'}
    if (wt_groups & primary_props) == (mut_groups & primary_props) and 'Proline' not in wt_groups and 'Methionine' not in wt_groups and 'Histidine' not in wt_groups:
        applicable_types.append(('Similar to WT', compound_to_initials(compound)))

    # Checks for Proline, Methionine and Histidine transitions
    if 'Proline' in wt_groups and 'Proline' not in mut_groups:
        applicable_types.append(('Proline-Any', compound_to_initials(compound)))
    if 'Proline' not in wt_groups and 'Proline' in mut_groups:
        applicable_types.append(('Any-Proline', compound_to_initials(compound)))
    if 'Methionine' in wt_groups and 'Methionine' not in mut_groups:
        applicable_types.append(('Methionine-Any', compound_to_initials(compound)))
    if 'Methionine' not in wt_groups and 'Methionine' in mut_groups:
        applicable_types.append(('Any-Methionine', compound_to_initials(compound)))
    if 'Histidine' in wt_groups and 'Histidine' not in mut_groups:
        applicable_types.append(('Histidine-Any', compound_to_initials(compound)))
    if 'Histidine' not in wt_groups and 'Histidine' in mut_groups:
        applicable_types.append(('Any-Histidine', compound_to_initials(compound)))


    # Check for WT transitions based on defined mutation transitions
    for transition, label in mutation_transitions.items():
        if (transition[0] in wt_groups and transition[1] in mut_groups):
            applicable_types.append((label, compound_to_initials(compound)))
    # Print for troubleshoot
    #print(f"Mutation details for position {position} ,WT aa {wt_aa} to Mut aa {mut_aa}: {applicable_types}")
    return applicable_types if applicable_types else [('Similar to WT', compound_to_initials(compound))]

# Applying mutation position/type functions
df_resistant['Position'] = df_resistant['aa_seq'].apply(find_mutation_position).astype(int)

# Applying mutation position and type functions
df_resistant['Mutation_Details'] = df_resistant.apply(
    lambda row: determine_mutation_type_and_resistance(
        wt_sequence, row['aa_seq'], row['Position'], row['compound']), axis=1)

# Print before the loop to inspect the details
#print("Mutation Details Example:", df_resistant['Mutation_Details'].iloc[0])

# Aggregation for pivot table
pivot_data = {}
for index, row in df_resistant.iterrows():
    position = row['Position']
    mutation_details = row['Mutation_Details']
    for mutation, initials in mutation_details:
        if mutation not in pivot_data:
            pivot_data[mutation] = {}
        if position not in pivot_data[mutation]:
            pivot_data[mutation][position] = set()
        pivot_data[mutation][position].update(initials.split(','))
df_resistant.to_csv('final_data_for_heatmap.csv', index=False)

# Convert nested dictionaries to dataframe for the heatmap
pivot_table = pd.DataFrame.from_dict({k: {pos: ''.join(sorted(vals)) for pos, vals in v.items()} for k, v in pivot_data.items()}, orient='index').fillna(' ')
pivot_table = pivot_table.reindex(ordered_mutation_types)  # Reorder the mutation types
pivot_table = pivot_table[sorted(pivot_table.columns)]  # Sort columns by residue position
pivot_table.fillna(' ', inplace=True)

# Function that checks if a mutation type is applicable based on WT amino acid groups
def is_transition_applicable(wt_groups, mutation_type, wt_aa):
    # Handle specific cases for Proline, Methionine and Histidine
    if mutation_type == 'Proline-Any' and 'Proline' not in wt_groups:
        return False
    if mutation_type == 'Any-Proline' and 'Proline' in wt_groups:
        return False
    if mutation_type == 'Methionine-Any' and 'Methionine' not in wt_groups:
        return False
    if mutation_type == 'Any-Methionine' and 'Methionine' in wt_groups:
        return False
    if mutation_type == 'Any-Proline' and 'Proline' not in wt_groups:
        return True
    if mutation_type == 'Any-Methionine' and 'Methionine' not in wt_groups:
        return True
    if mutation_type == 'Histidine-Any' and 'Histidine' not in wt_groups:
        return False
    if mutation_type == 'Any-Histidine' and 'Histidine' in wt_groups:
        return False
    if mutation_type == 'Any-Histidine' and 'Histidine' not in wt_groups:
        return True
    if mutation_type == 'Similar to WT':
        if 'Proline' in wt_groups or 'Methionine' in wt_groups or 'Histidine' in wt_groups:
            return False

    # General rule for transitions based on mutation properties
    parts = mutation_type.split('-')
    if len(parts) == 2:
        from_group, to_group = parts
        # Applicable if direct transition matches WT group properties
        if from_group in wt_groups and mutation_transitions.get((from_group, to_group), '') == mutation_type:
            return True
        # Not applicable if reverse transition matches WT group properties
        if to_group in wt_groups and mutation_transitions.get((to_group, from_group), '') == mutation_type:
            return False
    if mutation_type == 'Similar to WT':
        if 'Proline' not in wt_groups or 'Methionine' not in wt_groups or 'Histidine' not in wt_groups:
            return True
    return False

# Before visualization, mark non-applicable types
for column in sorted(pivot_table.columns):
    wt_aa = wt_sequence[int(column) - 1]  # Adjust indexing to match sequence positions, ensuring integer indexing
    wt_groups = set(group for group, aas in aa_properties.items() if wt_aa in aas)
    for mutation_type in pivot_table.index:
        if not is_transition_applicable(wt_groups, mutation_type, wt_aa):
            pivot_table.at[mutation_type, column] = '-'  # Mark as not applicable

# Printing to inspect the final pivot table
#print("Pivot Table after applicability check:")
#print(pivot_table)


# Visualization using matplotlib
fig, ax = plt.subplots(figsize=(5, 4))
cmap = ListedColormap(['whitesmoke', 'gold', 'tomato'])
data = pivot_table

# Draw each cell according to its status
for (i, j), val in np.ndenumerate(pivot_table.values):
    label = pivot_table.index[i]
    if val == '-':
        color = 'gainsboro'
        text = ''
    else:
        color = 'gold' if val.strip() else 'tomato'
        text = val
    ax.text(j + 0.5, i + 0.4, text, ha='center', va='center', color='black')
    ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=True, edgecolor='ivory', facecolor=color))

ax.set_xlim(0, len(pivot_table.columns))
ax.set_ylim(0, len(pivot_table.index))
ax.set_xticks(np.arange(len(pivot_table.columns)) + 0.5)
ax.set_xticklabels(list(wt_sequence))
ax.set_yticks(np.arange(len(pivot_table.index)) + 0.5)
ax.set_yticklabels(pivot_table.index, rotation=0)
ax.set_xlabel('WT Residues 639 to 647')
ax.set_ylabel('Mutation Type')
plt.grid(visible=False)
plt.tight_layout()
plt.savefig('hs1_resistance_heatmap.svg', format='svg', dpi=1000)
plt.show()
