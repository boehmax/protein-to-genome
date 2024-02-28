import subprocess # Used in download_genome_data and unzip_downloaded_files functions to run shell commands.
import os # Used in create_summary_file function to list files in a directory and join paths.
from Bio import Entrez # Used in retrieve_protein_info function to fetch IPG files using Entrez.
import pandas as pd #Used throughout your script for data manipulation and analysis (creating DataFrames, reading CSV files, concatenating DataFrames, etc.)
import numpy as np # Used in ipg_xml_to_dataframe function to create a range for iteration.

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


def generate_protein_alias_ipg():
    """
    From the IPG files, generate a file containing the protein alias.
    """
    # create a new file for the outpu   t
    output_file = "output/ipg_representative.txt"
    alias_df = pd.DataFrame()
    # loop over all .txt files in the directory
    for filename in os.listdir('ipg/*.csv'):
        # read the file into a pandas dataframe
        df = pd.read_csv(filename)

        # get the name of the first entry
        first_entry = df['protaccver'][0]

        # get the 'Protein' column and add the first entry as a new column
        df_output = df['protaccver'].copy()
        df_output['PIGI'] = first_entry

        # append the output dataframe to the overall dataframe
        alias_df = pd.concat([alias_df,df_output], ignore_index=True)

    # write the overall dataframe to the output file
    alias_df.to_csv(output_file, header=True, index=False)
    

def create_summary_file():
    """
    Create a summary file containing the second line of each IPG file.
    """
    summary_data = pd.DataFrame()
    for filename in os.listdir('ipg/'):
        with open(os.path.join('ipg/', filename), 'r') as file:
            reader = pd.read_csv(file)
            row1 = reader.head(1)
            summary_data = pd.concat([summary_data,row1], ignore_index=True)
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


