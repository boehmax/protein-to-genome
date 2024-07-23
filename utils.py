import subprocess # Used in download_genome_data and unzip_downloaded_files functions to run shell commands.
import os # Used in create_summary_file function to list files in a directory and join paths.
from Bio import Entrez # Used in retrieve_protein_info function to fetch IPG files using Entrez.
import pandas as pd #Used throughout your script for data manipulation and analysis (creating DataFrames, reading CSV files, concatenating DataFrames, etc.)
import numpy as np # Used in ipg_xml_to_dataframe function to create a range for iteration.
from io import StringIO



def ipg_xml_to_dataframe(ipg_xml):
    """
    Takes IPG output XML to make a dataframe of the results
    """
    ipg_record = pd.DataFrame()
    protein_list = ipg_xml.copy()['IPGReport']['ProteinList']
    # Iterate over each protein in the protein list
    for index in np.arange(len(protein_list)):
        buffer = protein_list[index]
        dict1 = buffer.attributes # getting the attributes of the first protein
        dict1['protaccver']=dict1.pop('accver', None)
        dict2 = buffer['CDSList'][0].attributes # getting the attributes of the first CDS
        dict2['nucaccver']=dict2.pop('accver', None)
        combined_dict = {**dict1, **dict2}
        # Check if data has a 'strain' key
        if 'strain' not in combined_dict:
            # If not, add the key with a None value
            combined_dict['strain'] = None
        combined_df = pd.DataFrame([combined_dict])
        ipg_record = pd.concat([ipg_record, combined_df], ignore_index=True)
    return ipg_record


def retrieve_protein_info(filename):
    """
    Retrieve IPG files for each protein ID listed in the input file.
    """
    with open(filename, 'r') as file:
        for line in file:
            protein_id = line.strip()
            output_file = f'ipg/{protein_id}.csv'
            # Fetch IPG file using Entrez
            stream = Entrez.efetch(db="protein", id=protein_id, rettype="ipg", retmode="xml")
            record = Entrez.read(stream)   # Read the response data
            ipg_data_df = ipg_xml_to_dataframe(record)
            with open(output_file, 'w') as ipg_file:
                ipg_data_df.to_csv(ipg_file, index=False, header=True)  # Write the decoded data to the file


def extend_ipg_files_with_assembly_information(directory='ipg'):
    """
    Extend IPG files with assembly length information.
    With NCBI datasets functionality.
    """
    for filename in os.listdir(directory):
        if filename.endswith('.csv'):
            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path)
            
            # Check if 'assembly' column exists
            if 'assembly' in df.columns:
                print("Processing file: ", file_path) 
                for index, row in df.iterrows():
                    assembly = row['assembly']
                    # Construct and execute the command
                    command = f"../../ncbi/datasets summary genome accession {assembly} --as-json-lines | ../../ncbi/dataformat tsv genome --fields accession,checkm-completeness,checkm-contamination,checkm-version,assmstats-contig-n50,assmstats-contig-l50,assmstats-total-ungapped-len,assmstats-total-sequence-len"
                    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                    stdout, stderr = process.communicate()
                    if process.returncode == 0:
                        # Assuming the output is in TSV format as specified
                        fields = pd.read_csv(StringIO(stdout), sep='\t')
                        if len(fields.columns) == 8:
                            # Update the DataFrame with the new information
                            df.at[index, 'assembly_accession'] = fields.iloc[0, 0]
                            df.at[index, 'checkm_completeness'] = fields.iloc[0, 1]
                            df.at[index, 'checkm_contamination'] = fields.iloc[0, 2]
                            df.at[index, 'checkm_version'] = fields.iloc[0, 3]
                            df.at[index, 'contig_N50'] = fields.iloc[0, 4]
                            df.at[index, 'contig_L50'] = fields.iloc[0, 5]
                            df.at[index, 'ungaped_seq_len'] = fields.iloc[0, 6]
                            df.at[index, 'seq_len'] = fields.iloc[0, 7]
                            print(df)
                    else:
                        print(f"Error processing {assembly}: {stderr}")
                #Sort the dataframe by source (RefSeq) and checkm completeness, contig N50, contig L50, ungaped sequence length, sequence length
                df = df.sort_values(by=['source', 'checkm_completeness', 'contig_N50', 'contig_L50', 'ungaped_seq_len', 'seq_len'], ascending=[False, False, False, True, False, False])
                
                # Write the updated DataFrame back to the file
                df.to_csv(file_path, index=False)
            else:
                print(f"'assembly' column not found in {file_path}")

def generate_protein_alias_ipg():
    """
    From the IPG files, generate a file containing the protein alias.
    """
    # create a new file for the outpu   t
    output_file = "output/ipg_representative.txt"
    alias_df = pd.DataFrame()
    directory = 'ipg'  # Directory containing the files

    # List all files in the directory and filter for .csv files
    for filename in os.listdir(directory):
        if filename.endswith('.csv'):
            # Construct the full path of the file
            file_path = os.path.join(directory, filename)
            # Read the file into a pandas dataframe
            df = pd.read_csv(file_path)

            # Get the name of the first entry
            first_entry = df['protaccver'][0]

            # Create a new DataFrame from 'protaccver' column
            df_output = pd.DataFrame(df['protaccver'])
            df_output['PIGI'] = first_entry  # Add the first entry as a new column

            # Append the output dataframe to the overall dataframe
            alias_df = pd.concat([alias_df, df_output], ignore_index=True)

    # Write the overall dataframe to the output file
    alias_df.to_csv(output_file, header=True, index=False)
    

def create_summary_file():
    """
    Create a summary file containing the second line of each IPG file.
    """
    summary_data = pd.DataFrame()
    for filename in os.listdir('ipg/'):
        if filename.endswith('.csv'):  # Check if the file is a .csv file
            with open(os.path.join('ipg/', filename), 'r') as file:
                reader = pd.read_csv(file)
                if len(reader) > 1:  # Ensure there is at least a second line
                    row1 = reader.iloc[1:2]  # Get the second line
                    summary_data = pd.concat([summary_data, row1], ignore_index=True)
    summary_data.to_csv('ipg_summary.csv', index=False, header=True)
 
def process_summary_file():
    """
    Process the summary file to extract relevant information.
    """
    summary_df = pd.read_csv('ipg_summary.csv') 
    accs_protein = summary_df[["assembly", "protaccver"]]
    accs = accs_protein['assembly']
    accs.to_csv('assm_accs.csv', index=False, header=False)


def download_genome_data():
    """
    Download genome data based on the list of assembly accessions.
    """
    PATH_TO_NCBI_DATASETS = '../../ncbi/datasets'
    subprocess.run([PATH_TO_NCBI_DATASETS, 'download', 'genome', 'accession', '--inputfile', 'assm_accs.csv', '--include', 'gff3'])

def unzip_downloaded_files():
    """
    Unzip the downloaded files.
    """
    subprocess.run(['unzip', 'ncbi_dataset.zip'])


