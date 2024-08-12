import pandas as pd

# Data for nucleotides
nucleotide_data = {
    'Nucleotide': ['DNA A (bz)', 'DNA C (bz)', 'DNA G (ibu)', 'DNA T', 'RNA A (bz)', 'RNA C (bz)', 'RNA G (ibu)', 'RNA U',
                   '2\'OMe A (bz)', '2\'OMe C (bz)', '2\'OMe G (ibu)', '2\'OMe U', '2\'F-ANA A', '2\'F-ANA C', 
                   '2\'F-ANA G', '2\'F-ANA T', '2\'F-ANA U', '2\'F A', '2\'F C', '2\'F G', '2\'F U', 
                   'Ara A', 'Ara C', 'Ara G', 'Ara U', 'LNA A', 'LNA C', 'LNA G', 'LNA T', '~DBCO-Serinol'],
    'Acronym': ['dA', 'dC', 'dG', 'dT', 'rA', 'rC', 'rG', 'rU', 'oA', 'oC', 'oG', 'oU', 'faA', 'faC', 'faG', 
                'faT', 'faU', 'fA', 'fC', 'fG', 'fU', 'aA', 'aC', 'aG', 'aU', 'lA', 'lC', 'lG', 'lT', 'DBCO'],
    'Molar Mass (g/mol)': [857.95, 833.93, 839.94, 744.83, 988.18, 964.17, 970.16, 861.03, 
                           887.958, 863.934, 869.941, 760.812, 875.92, 851.9, 857.91, 762.803, 748.78, 875.92, 
                           851.9, 857.91, 748.78, 915.97, 891.94, 897.95, 788.82, 885.94, 861.92, 852.91, 909.08, 909.08],
    'Exact Mass': [313.058, 289.046, 329.053, 304.046, 329.053, 305.041, 345.047, 306.025, 
                   343.07, 319.06, 359.06, 320.04, 331.0486, 307.0374, 347.0435, 322.0371, 308.0371, 331.0486, 
                   307.0374, 347.0435, 308.0371, 313.058, 289.046, 329.053, 304.046, 341.053, 317.0375, 357.048, 332.0415, 406.4912],
    'Avg Mass': [313.207, 289.182, 329.206, 304.193, 329.206, 305.181, 345.205, 306.166, 
                 343.24, 319.21, 359.24, 320.19, 331.1975, 307.1728, 347.1969, 322.1842, 308.258, 331.1975, 
                 307.1728, 347.1969, 308.258, 313.207, 289.182, 329.206, 304.193, 341.2171, 317.199, 357.2165, 332.2038, 406.4853],
    'Coupling Time (s)': [180, 180, 240, 180, 720, 720, 780, 720, 360, 360, 420, 360, 180, 180, 240, 180, 180, 360, 360, 420, 360, 720, 720, 780, 720, 180, 180, 240, 180, 720]
}

# Create a DataFrame
df = pd.DataFrame(nucleotide_data)

# Sequences
sequences = {
    "Order": "DBCO faU faC dA dC dG faT"
}

# Calculation constants
uL_per_insert = 300
molarity = 0.1
extra_inserts = 5

# Prepare results lists
results = []
sequence_results = []

# Process each sequence
for seq_name, sequence in sequences.items():
    nucleotides = sequence.split()

    # Count occurrences and maintain order
    nucleotide_counts = {}
    exact_mass = 0
    avg_mass = 0
    ps_linkage_count = 0

    for nucleotide in nucleotides:
        base_acronym = nucleotide.strip('*')
        if base_acronym in nucleotide_counts:
            nucleotide_counts[base_acronym] += 1
        else:
            nucleotide_counts[base_acronym] = 1

        matched_rows = df[df['Acronym'] == base_acronym]
        if not matched_rows.empty:
            matched_row = matched_rows.iloc[0]
            exact_mass += matched_row['Exact Mass']
            avg_mass += matched_row['Avg Mass']
            ps_linkage_count += nucleotide.count('*')
        else:
            print(f"Warning: Acronym '{base_acronym}' not found in nucleotide data")

    # Add PS linkage adjustments
    exact_mass += ps_linkage_count * 16
    avg_mass += ps_linkage_count * 16

    # Subtract 5' OH adjustment
    exact_mass -= 61.9558
    avg_mass -= 61.9647

    # Process each nucleotide to handle the dash and star properly
    formatted_nucleotides = []
    for nucleotide in nucleotides:
        if '*' in nucleotide:
            formatted_nucleotides.append(nucleotide.replace('*', '^'))
        else:
            formatted_nucleotides.append(nucleotide)
    
    # Join the nucleotides with a dash, ensuring correct formatting
    formatted_sequence = "5'-" + '-'.join(formatted_nucleotides).replace('^-', '^').replace('^-', '^') + "-3'"

    # Create sequence results
    sequence_results.append([
        seq_name, formatted_sequence, len(nucleotides), exact_mass, avg_mass
    ])

    # Perform calculations for each nucleotide
    for base_acronym, count in nucleotide_counts.items():
        matched_rows = df[df['Acronym'] == base_acronym]
        if not matched_rows.empty:
            matched_row = matched_rows.iloc[0]
            total_inserts = count + extra_inserts
            mg = matched_row['Molar Mass (g/mol)'] * molarity * (uL_per_insert * total_inserts / 1000)
            mL = (uL_per_insert * total_inserts / 1000)
            
            results.append([
                matched_row['Nucleotide'], base_acronym, matched_row['Molar Mass (g/mol)'],
                matched_row['Coupling Time (s)'],
                count, extra_inserts, total_inserts, uL_per_insert, molarity, mg, mL
            ])

# Combine results for the same nucleotide
results_df = pd.DataFrame(results, columns=[
    'Nucleotide', 'Acronym', 'Molar Mass (g/mol)', 'Coupling Time (s)', 'Sequence Inserts', 'Extra Inserts', 
    'Total Inserts', 'uL per Insert', 'M (Molarity)', 'mg', 'mL'
])

combined_results_df = results_df.groupby(['Nucleotide', 'Acronym', 'Molar Mass (g/mol)', 'Coupling Time (s)']).agg(
    {'Sequence Inserts': 'sum', 'Extra Inserts': 'max', 'Total Inserts': 'sum', 'uL per Insert': 'max', 'M (Molarity)': 'max', 'mg': 'sum', 'mL': 'sum'}
).reset_index()

# Calculate Total Inserts using the aggregated values
combined_results_df['Total Inserts'] = combined_results_df['Sequence Inserts'] + combined_results_df['Extra Inserts']

sequence_results_df = pd.DataFrame(sequence_results, columns=[
    'Sequence Name', 'Sequence', 'Oligo Length', 'Exact Mass', 'Average Mass'
])

# Function to modify the Mermade Code based on asterisk
def modify_mermade_code(nucleotide, mermade_char):
    if '*' in nucleotide:
        if mermade_char.isdigit():
            return str(int(mermade_char) - 4)
        else:
            return mermade_char.lower()
    return mermade_char

# Function to generate the Mermade codes for each sequence
def generate_mermade_codes(sequences, acronym_to_mermade):
    sequence_codes = []
    for seq_name, sequence in sequences.items():
        nucleotides = sequence.split()
        formatted_sequence = sequence_results_df[sequence_results_df['Sequence Name'] == seq_name]['Sequence'].values[0]
        mermade_code = ''.join(modify_mermade_code(n, acronym_to_mermade[n.strip('*')]) for n in nucleotides)
        mermade_code_formatted = f"{seq_name}, {mermade_code}DBB"
        sequence_codes.append([seq_name, formatted_sequence, mermade_code_formatted])
    return pd.DataFrame(sequence_codes, columns=['Sequence Name', 'Sequence', 'Mermade Code'])

# Define the mapping for Mermade Code
mermade_mapping = {
    1: 'A', 2: 'C', 3: 'G', 4: 'U', 5: '5', 6: '6', 7: '7', 8: '8', 9: 'W', 10: 'X', 11: 'Y', 12: 'Z', 13: 'NA', 14: 'NA', 15: 'NA', 16: 'NA'
}

# Reset index due to 'Acronym' error
combined_results_df = results_df.groupby(['Nucleotide', 'Acronym', 'Molar Mass (g/mol)', 'Coupling Time (s)']).agg(
    {'Sequence Inserts': 'sum', 'Extra Inserts': 'max', 'Total Inserts': 'sum', 
     'uL per Insert': 'max', 'M (Molarity)': 'max', 'mg': 'sum', 'mL': 'sum'}
).reset_index()

# Create a dictionary to map Acronym to Mermade Code
acronym_to_mermade = {acronym: mermade_mapping[idx+1] for idx, acronym in enumerate(combined_results_df['Acronym'].unique())}
unique_acronyms = combined_results_df['Acronym'].unique()

# Add Mermade Code column based on the line number
combined_results_df.insert(3, 'Mermade Code', combined_results_df['Acronym'].map(acronym_to_mermade))

# Generate the Mermade codes for each sequence
mermade_codes_df = generate_mermade_codes(sequences, acronym_to_mermade)

# Function to convert sequences to IDT Code
def convert_to_idt_code(sequence):
    sequence = sequence.replace("5'-", "").replace("-3'", "")
    idt_code = []

    nucleotides = sequence.split('-')
    for i, nucleotide in enumerate(nucleotides):
        if nucleotide == 'DBCO':
            idt_code.append('/5DBCON/')
        elif i == 0:  # 5' terminal modification
            if nucleotide.startswith('fa'):
                idt_code.append('/52F' + nucleotide[2] + '/')
            elif nucleotide.startswith('f'):
                idt_code.append('/52F' + nucleotide[1] + '/')
        elif i == len(nucleotides) - 1:  # 3' terminal modification
            if nucleotide.startswith('fa'):
                idt_code.append('/32F' + nucleotide[2] + '/')
            elif nucleotide.startswith('f'):
                idt_code.append('/32F' + nucleotide[1] + '/')
        else:
            if nucleotide.startswith('fa'):
                idt_code.append('/i2F' + nucleotide[2] + '/')
            elif nucleotide.startswith('f'):
                idt_code.append('/i2F' + nucleotide[1] + '/')
            elif nucleotide.startswith('d'):
                idt_code.append(nucleotide[1:])
            elif nucleotide.startswith('o'):
                idt_code.append('m' + nucleotide[1:])
            elif nucleotide.startswith('l'):
                idt_code.append('+' + nucleotide[1:])
            else:
                idt_code.append(nucleotide)
    
    return ''.join(idt_code)

# Function to generate the IDT codes for each sequence
def generate_idt_codes(sequences):
    sequence_codes = []
    for seq_name, sequence in sequences.items():
        idt_code = convert_to_idt_code("5'-" + '-'.join(sequence.split()) + "-3'")
        sequence_codes.append([seq_name, sequence, idt_code])
    return pd.DataFrame(sequence_codes, columns=['Sequence Name', 'Sequence', 'IDT Code'])

# Generate the IDT codes for each sequence
idt_codes_df = generate_idt_codes(sequences)

# Merge Mermade and IDT codes
sequence_codes_df = pd.merge(mermade_codes_df, idt_codes_df[['Sequence Name', 'IDT Code']], on='Sequence Name')

# Alternative protecting groups data
alternative_data = {
    'Nucleotide': ['dC (ac)', 'dG (ac)', 'rC (ac)', 'rG (ac)', '2-OMe-rG (PAC)'],
    'Acronym': ['dC', 'dG', 'rC', 'rG', 'oG'],
    'Molar Mass (g/mol)': [771.84, 811.86, 902.1, 942.12, 933.96]
}

# Create a DataFrame for the alternative protecting groups
alternative_df = pd.DataFrame(alternative_data)
    
    # Function to create alternative protecting groups DataFrame with calculations
def create_alternative_protecting_groups_sheet(main_df, alternative_df, uL_per_insert, molarity, extra_inserts):
    results = []
    
    for _, alt_row in alternative_df.iterrows():
        acronym = alt_row['Acronym']
        matched_rows = main_df[main_df['Acronym'] == acronym]
        
        if not matched_rows.empty:
            for _, matched_row in matched_rows.iterrows():
                total_inserts = matched_row['Sequence Inserts'] + extra_inserts
                mL = (uL_per_insert * total_inserts) / 1000
                mg = alt_row['Molar Mass (g/mol)'] * molarity * mL
                
                results.append([
                    alt_row['Nucleotide'], acronym, alt_row['Molar Mass (g/mol)'], 
                    matched_row['Coupling Time (s)'], 
                    matched_row['Mermade Code'], matched_row['Sequence Inserts'], 
                    extra_inserts, total_inserts, uL_per_insert, molarity, mg, mL
                ])
    
    return pd.DataFrame(results, columns=[
        'Nucleotide', 'Acronym', 'Molar Mass (g/mol)', 'Coupling Time (s)', 'Mermade Code', 
        'Sequence Inserts', 'Extra Inserts', 'Total Inserts', 
        'uL per Insert', 'M (Molarity)', 'mg', 'mL'
    ])

# Calculate alternative protecting groups
alternative_protecting_groups_df = create_alternative_protecting_groups_sheet(combined_results_df, alternative_df, uL_per_insert, molarity, extra_inserts)

# Save the DataFrames to an Excel file with auto-adjusted column widths
with pd.ExcelWriter('oligo_calculations.xlsx', engine='openpyxl') as writer:
    combined_results_df.to_excel(writer, sheet_name='Nucleotide Table', index=False)
    alternative_protecting_groups_df.to_excel(writer, sheet_name='Alternative Protecting Groups', index=False)
    sequence_results_df.to_excel(writer, sheet_name='Sequence Table', index=False)
    sequence_codes_df.to_excel(writer, sheet_name='Sequence Codes', index=False)
    
    workbook = writer.book
    for sheetname in writer.sheets:
        worksheet = writer.sheets[sheetname]
        for col in worksheet.columns:
            max_length = 0
            column = col[0].column_letter # Get the column name
            for cell in col:
                try:
                    if len(str(cell.value)) > max_length:
                        max_length = len(cell.value)
                except:
                    pass
            adjusted_width = (max_length + 2)
            worksheet.column_dimensions[column].width = adjusted_width
            
    # Adjust the width for 'mg' and 'mL' columns
    worksheet = writer.sheets['Nucleotide Table']
    worksheet.column_dimensions['K'].width = 20  # 'mg' column
    worksheet.column_dimensions['L'].width = 20  # 'mL' column
    worksheet = writer.sheets['Alternative Protecting Groups']
    worksheet.column_dimensions['K'].width = 20  # 'mg' column
    worksheet.column_dimensions['L'].width = 20  # 'mL' column

print("Excel file 'Oligo Calculations.xlsx' has been created with the combined tables.")