Welcome!

This repository contains a Python script that performs oligonucleotide synthesis calculations, generates IDT and Mermade codes, and exports the results to an Excel file. 
The script handles various modifications in nucleotide sequences and ensures proper formatting for oligonucleotide synthesis.

Features:
Nucleotide Table Generation: Creates a table of nucleotides with their molar mass, exact mass, average mass, and coupling times.
Sequence Processing: Processes given nucleotide sequences to calculate exact and average mass, taking into account any phosphorothioate linkages.
IDT Code Generation: Converts sequences into IDT-compatible codes, applying specific rules for nucleotide modifications, including DBCO conversion and terminal modifications.
Mermade Code Generation: Generates Mermade codes for sequences based on custom mappings.
Alternative Protecting Groups: Provides calculations for alternative protecting groups.
Excel Output: Outputs the results to an Excel file, applying borders to all cells and adjusting column widths for readability.

Running the Script:
Place your sequences in the sequences dictionary, with each acronym separated by a space. For PS linkages, add a '*' after the acronym.
Run the script. The script will process each sequence, perform necessary calculations, and generate IDT and Mermade codes.
The results are saved in an Excel file named Oligo Calculations.xlsx.

Best of luck!
