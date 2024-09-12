import subprocess # Used in download_genome_data and unzip_downloaded_files functions to run shell commands.
import os # Used in create_summary_file function to list files in a directory and join paths.
from Bio import Entrez # Used in retrieve_protein_info function to fetch IPG files using Entrez.
import pandas as pd #Used throughout your script for data manipulation and analysis (creating DataFrames, reading CSV files, concatenating DataFrames, etc.)
import numpy as np # Used in ipg_xml_to_dataframe function to create a range for iteration.
from io import StringIO
from datetime import datetime

# Compute the current date once
current_date = datetime.now().strftime('%Y-%m-%d')

def get_current_date_directory(base_directory='output', subdirectory='ipg', date = current_date) -> str:
    """
    Get the current date directory in the format YYYY-MM-DD.
    
    Parameters:
    - base_directory (str): The base directory to store the files.
    - subdirectory (str): The subdirectory to store the files.
    - date (str): The current date in the format YYYY-MM-DD.
    
    Returns:
    - directory (str): The full path to the current date
    """
    current_date = date
    directory = os.path.join(base_directory, current_date, subdirectory)
    os.makedirs(directory, exist_ok=True)  # Create the directory if it doesn't exist
    return directory

# Function to safely extract attributes from StringElement needed for XML parsing
def extract_attributes(element)  -> dict:
    """
    Safely extract attributes from an XML element.

    Parameters:
    - element: The XML element.

    Returns:
    - dict: The attributes of the element.
    """
    if hasattr(element, 'attributes'):
        return element.attributes
    return {}


def merge_dicts_with_suffix(*dicts, suffix="_dup"):
    merged_dict = {}
    for d in dicts:
        for key, value in d.items():
            new_key = key
            while new_key in merged_dict:
                new_key += suffix
            merged_dict[new_key] = value
    return merged_dict

# Function to convert the dictionary into a DataFrame needed for XML parsing
def dict_to_df(data) -> pd.DataFrame:
    """
    Convert a dictionary into a DataFrame.

    Parameters:
    - data (dict): The dictionary to convert.

    Returns:
    - pd.DataFrame: The resulting DataFrame.
    """
    ipg_report = data.get('IPGReport', {})
    
    # Extract Product Information
    product_info = extract_attributes(ipg_report.get('Product', {}))
    
    # Extract ProteinList Information
    protein_list = ipg_report.get('ProteinList', [])
    protein_data = []
    for protein in protein_list:
        protein_attrs = extract_attributes(protein)
        cds_list = protein.get('CDSList', [])
        for cds in cds_list:
            row = merge_dicts_with_suffix(product_info, protein_attrs, extract_attributes(cds))
            protein_data.append(row)
    
    # Convert to DataFrame
    df = pd.DataFrame(protein_data)
    
    return df
    

def retrieve_protein_info(filename) -> None:
    """
    Retrieve IPG files for each protein ID listed in the input file.

    Parameters:
    - filename (str): The path to the input file containing protein IDs.
    """
    directory = get_current_date_directory()
    print("Retrieving IPG files...")
    with open(filename, 'r') as file:
        for line in file:
            protein_id = line.strip()
            output_file = f'{directory}/{protein_id}.csv'
            # Fetch IPG file using Entrez
            stream = Entrez.efetch(db="protein", id=protein_id, rettype="ipg", retmode="xml")
            record = Entrez.read(stream)   # Read the response data
            ipg_data_df = dict_to_df(record)
            # Add PIGI column
            ipg_data_df['PIGI'] = protein_id
            print(f"Retrieved IPG data for {protein_id}.")
            with open(output_file, 'w') as ipg_file:
                ipg_data_df.to_csv(ipg_file, index=False, header=True)  # Write the decoded data to the file
    print("IPG files retrieved.")


def extend_ipg_files_with_assembly_information() -> None:
    """
    Extend IPG files with assembly length information.
    With NCBI datasets functionality.
    """
    directory = get_current_date_directory()
    print("Extending IPG files with assembly information...")
    for filename in os.listdir(directory):
        if filename.endswith('.csv'):
            file_path = os.path.join(directory, filename)
            
            # Check if the file is empty before reading
            if os.path.getsize(file_path) > 0:
                try:
                    df = pd.read_csv(file_path)
                except pd.errors.EmptyDataError:
                    print(f"Error processing {file_path}: No columns to parse from file.")
                    continue  # Skip to the next file
            else:
                print(f"Error processing {file_path}: File is empty.")
                continue  # Skip to the next file
            
            # Check if 'assembly' column exists
            if 'assembly' in df.columns:
                print("Processing file: ", file_path) 
                for index, row in df.iterrows():
                    assembly = row['assembly']
                    # Construct and execute the command
                    command = f"~/ncbi/datasets summary genome accession {assembly} --as-json-lines | ~/ncbi/dataformat tsv genome --fields accession,checkm-completeness,checkm-contamination,checkm-version,assmstats-contig-n50,assmstats-contig-l50,assmstats-total-ungapped-len,assmstats-total-sequence-len,source_database"
                    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                    stdout, stderr = process.communicate()
                    if process.returncode == 0:
                        if stdout.strip():  # Check if stdout is not empty or just whitespace
                            try:
                                fields = pd.read_csv(StringIO(stdout), sep='\t')
                            except pd.errors.EmptyDataError:
                                print(f"Error processing {assembly}: No columns to parse from stdout.")
                                continue  # Skip to the next row

                            if len(fields.columns) == 9:
                                # Update the DataFrame with the new information
                                # List of required columns
                                required_columns = [
                                    'assembly_accession', 'checkm_completeness', 'checkm_contamination', 
                                    'checkm_version', 'contig_N50', 'contig_L50', 'ungaped_seq_len', 'seq_len', 'source'
                                ]

                                # Ensure all required columns exist in the DataFrame
                                for column in required_columns:
                                    if column not in df.columns:
                                        df[column] = None

                                # Assign values to the DataFrame        
                                df.at[index, 'assembly_accession'] = fields.iloc[0, 0]
                                df.at[index, 'checkm_completeness'] = fields.iloc[0, 1]
                                df.at[index, 'checkm_contamination'] = fields.iloc[0, 2]
                                df.at[index, 'checkm_version'] = fields.iloc[0, 3]
                                df.at[index, 'contig_N50'] = fields.iloc[0, 4]
                                df.at[index, 'contig_L50'] = fields.iloc[0, 5]
                                df.at[index, 'ungaped_seq_len'] = fields.iloc[0, 6]
                                df.at[index, 'seq_len'] = fields.iloc[0, 7]
                                df.at[index, 'source'] = fields.iloc[0, 8]
                            else:
                                print(f"Error processing {assembly}: Unexpected number of fields ({len(fields.columns)}) in the output.")
                        else:
                            print(f"Error processing {assembly}: No data in stdout.")
                    else:
                        print(f"Error processing {assembly}: {stderr}")

                # Sort the dataframe by source (RefSeq) and checkm completeness, contig N50, contig L50, ungaped sequence length, sequence length
                try:
                    df = df.sort_values(by=['source', 'checkm_completeness', 'contig_N50', 'contig_L50', 'ungaped_seq_len', 'seq_len'], ascending=[False, False, False, True, False, False])
                except KeyError as e:
                    print(f"KeyError: {e} - One or more columns are missing from the DataFrame.")
                
                # Write the updated DataFrame back to the file
                df.to_csv(file_path, index=False)
            else:
                print(f"'assembly' column not found in {file_path}")
    print("IPG files extended with assembly information.")

def generate_protein_alias_ipg() -> None:
    """
    From the IPG files, generate a file containing the protein alias.
    """
    output_directory = get_current_date_directory(subdirectory='summary')
    # Create a new file for the output
    output_file = f"{output_directory}/ipg_representative.txt"
    alias_df = pd.DataFrame()
    directory = get_current_date_directory()  # Directory containing the files

    # List all files in the directory and filter for .csv files
    for filename in os.listdir(directory):
        if filename.endswith('.csv'):
            # Construct the full path of the file
            file_path = os.path.join(directory, filename)
            filename_without_extension = os.path.splitext(os.path.basename(file_path))[0]
            
            # Check if the file is empty before reading
            if os.path.getsize(file_path) > 0:
                try:
                    df = pd.read_csv(file_path)
                except pd.errors.EmptyDataError:
                    print(f"Error processing {file_path}: No columns to parse from file.")
                    continue  # Skip to the next file
            else:
                print(f"Error processing {file_path}: File is empty.")
                continue  # Skip to the next file

            # Check if 'accver' column exists
            if 'accver_dup' in df.columns:
                
                # Create a new DataFrame from 'accver' column
                df_output = pd.DataFrame(df['accver_dup'])
                df_output['PIGI'] = filename_without_extension  # Add the name of the file as a column

                # Append the output dataframe to the overall dataframe
                alias_df = pd.concat([alias_df, df_output], ignore_index=True)
            else:
                print(f"'accver_dup' column not found in {file_path}")

    # Write the overall dataframe to the output file
    alias_df.to_csv(output_file, header=True, index=False)

def create_summary_file() -> None:
    """
    Create a summary file containing the second line of each IPG file.
    """
    summary_data = pd.DataFrame()
    output_directory = get_current_date_directory(subdirectory='summary')
    directory = get_current_date_directory(subdirectory='ipg')  # Directory containing the files

    # List all files in the directory and filter for .csv files
    for filename in os.listdir(directory):
        if filename.endswith('.csv'):
            # Construct the full path of the file
            file_path = os.path.join(directory, filename)
            
            # Check if the file is empty before reading
            if os.path.getsize(file_path) > 0:
                try:
                    reader = pd.read_csv(file_path)
                except pd.errors.EmptyDataError:
                    print(f"Error processing {file_path}: No columns to parse from file.")
                    continue  # Skip to the next file
            else:
                print(f"Error processing {file_path}: File is empty.")
                continue  # Skip to the next file

            # Ensure there is at least a second line
            if len(reader) > 0:
                row1 = reader.iloc[[0]].copy()  # Select the first row and ensure it's a DataFrame
                row1.loc[:, 'PIGI'] = filename[:-4]  # Add the PIGI column and fill with file name without .csv  
                summary_data = pd.concat([summary_data, row1], ignore_index=True)
                #Add the PIGI column and ill with file name without .csv
                
            else:
                print(f"File {file_path} does not have any line.")

    # Write the overall dataframe to the output file
    summary_data.to_csv(f'{output_directory}/ipg_summary.csv', index=False, header=True)
 
def process_summary_file() -> None:
    """
    Process the summary file to extract relevant information.
    """
    output_directory = get_current_date_directory(subdirectory='summary')
    summary_df = pd.read_csv(f'{output_directory}/ipg_summary.csv') 
    accs_protein = summary_df[["assembly", "PIGI"]]
    accs = accs_protein['assembly']
    # Replace empty strings with NaN
    accs_protein = accs_protein.replace('', np.nan)
    accs = accs.replace('', np.nan)
    # Drop rows with NaN values to ensure only valid assembly accessions are used
    accs_protein = accs_protein.dropna()
    accs = accs.dropna()
    # Write the assembly accessions to a file
    accs_protein.to_csv(f'{output_directory}/assm_accs_protein.csv', index=False, header=False)
    accs.to_csv(f'{output_directory}/assm_accs.csv', index=False, header=False)

def download_genome_data() -> None:
    """
    Download genome data based on the list of assembly accessions.
    """
    path = os.getcwd()
    output_direcotry = get_current_date_directory(subdirectory='summary')
    os.chdir(output_direcotry)
    print("Downloading genome data...")
    PATH_TO_NCBI_DATASETS = os.path.expanduser('~/ncbi/datasets')
    subprocess.run([PATH_TO_NCBI_DATASETS, 'download', 'genome', 'accession', '--inputfile', 'assm_accs.csv', '--include', 'gff3'])
    os.chdir(path)

def unzip_downloaded_files() -> None:
    """
    Unzip the downloaded files.
    """
    path = os.getcwd()
    output_directory = get_current_date_directory(subdirectory='summary')
    os.chdir(output_directory)
    print("Extracting downloaded files...")
    subprocess.run(['unzip', 'ncbi_dataset.zip'])
    os.chdir(path)
    

