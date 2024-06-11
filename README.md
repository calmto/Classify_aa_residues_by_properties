# Classify_aa_residues_by_properties
A Python3 script from our study on protein-drug interactions. We analyze single-point mutations in a protein's drug-binding segments, categorizing mutations by amino acid properties ('Aromatic', 'Aliphatic', 'Positive', etc.). The resistance profiles are visualized in a heatmap using matplotlib, aiding in understanding drugs sensitivity mechanisms.

________________________________________
Title: Analysis of a protein's drug binding site's tolerance/sensitivity by deep mutational scanning (DMS) and classification of amino acid properties into 'mutation types'
________________________________________
Description: This script requires a hardcoded wild-type (WT) sequence and a CSV file in the current directory. The CSV file should contain the DMS results combined with the drug resistance screening of all mutants with the format: ,compound,seq_type,Nham_aa,aa_seq,s,cscore,refined_class,sensres 1,anidulafungin,single,0.0,FLVLSLRDP,0.06276288377908285,1.0,WT-like,sensitive 2,anidulafungin,single,1.0,*LVLSLRDP,-0.2306350142092216,1.0,WT-like,sensitive 3,anidulafungin,single,1.0,ALVLSLRDP,1.7562999789995668,1.0,intermediary,resistant 4,caspofungin,single,0.0,FLVLSLRDP,0.08276288377908285,1.0,WT-like,sensitive
This script was developed as part of a research project on preventive drug design of drug-resistant pathogens. The analysis involves: 1) creating mutant strains of every residue in the protein's drug binding site through DMS, 2) measure how the mutant strains react to the drug(s) compared to a WT strain (degree of resistance/tolerance to a drug), 3) Classify resistant mutations based on the amino acid properties of the WT and mutant residue, 4) aggregate and plot the dataframe to a heatmap for easy visualization of drug-specific resistance mechanism.
Useful for analyzing large datasets of single point mutations in a WT sequence based on a measured feature (resistance/tolerance in this case) and aa properties.
For more details, refer to the corresponding scientific article: [Link to the article]
________________________________________
Written by: Alexandre Torbey, PhD Candidate, Department of Health and Biotechnology, Armand Frappier Institute, Quebec (Canada). Main supervisors: Prof. David Chatenet (INRS-IAF) and Prof. Patrick Lague (Laval University) Co-authors: Dr. Romain Durand and Prof. Christian Landry (Laval University) Date: 11/June/2024
________________________________________
Dependencies:
•	numpy
•	pandas
•	matplotlib
To install the necessary dependencies, run the following command: pip install numpy pandas matplotlib
Ran on Python 3.10.10
________________________________________
License: This project is licensed under the MIT License - see the LICENSE file for details.

