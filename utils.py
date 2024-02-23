import subprocess
import os
from Bio import Entrez
import pandas as pd
from io import StringIO

def ipg_xml_to_dataframe(ipg_xml):
    """
    Takes IPG output XML to make a dataframe of the results
    """
    protein_list = ipg_xml['IPGReport']['ProteinList']
    cds_data = []
    # Iterate over each protein in the protein list
    for protein in protein_list:
        # Each protein has a CDSList which is a list of StringElement objects
        # We convert each StringElement to a dictionary and add it to our list
        for string_element in protein['attriubutes']:
            # Convert the StringElement to a dictionary
            dict_element = vars(string_element)
            # Pop out atributes filed
            #protein_attributes = dict_element.pop('attributes')
            # Add the dictionary to our list
            #cds_data.append(protein_attributes)
            cds_data.append(dict_element)

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(cds_data)

    return df

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
                ipg_data_df.to_csv(output_file, index=False, header=False)  # Write the decoded data to the file


def create_summary_file():
    """
    Create a summary file containing the second line of each IPG file.
    """
    summary_data = []
    for filename in os.listdir('ipg/'):
        with open(os.path.join('ipg/', filename), 'r') as file:
            reader = pd.read_csv(f)
            row1 = reader[0] 
            summary_data.append(row1)
    summary_df = pd.DataFrame(summary_data)
    summary_df.columns = ['Id',	'Source', 'Nucleotide Accession', 'Start', 'Stop', 'Strand', 'Protein', 'Protein Name','Organism','Strain','Assembly']
    summary_df.to_csv('ipg_summary.csv', index=False, header=False)
 
def process_summary_file():
    """
    Process the summary file to extract relevant information.
    """
    summary_df = pd.read_csv('ipg_summary.csv') 
    accs_protein = summary_df[["Assembly", "Protein"]]
    accs_protein.columns = ['Accession', 'Protein']
    accs_protein.to_csv('assm_accs_protein.csv', index=False)

    accs = accs_protein['Accession']
    accs[accs.str.startswith('GC')].to_csv('assm_accs.csv', index=False)


def download_genome_data():
    """
    Download genome data based on the list of assembly accessions.
    """
    PATH_TO_NCBI_DATASETS = '../../ncbi/datasets'
    subprocess.run([PATH_TO_NCBI_DATASETS, 'download', 'genome', 'accession', '--inputfile', 'assm_accs.txt', '--include', 'gff3'])

def unzip_downloaded_files():
    """
    Unzip the downloaded files.
    """
    subprocess.run(['unzip', 'ncbi_dataset.zip'])


