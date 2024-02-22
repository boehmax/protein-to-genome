import subprocess
import os
from Bio import Entrez
import pandas as pd
from io import StringIO

def retrieve_protein_info(filename):
    """
    Retrieve IPG files for each protein ID listed in the input file.
    """
    with open(filename, 'r') as file:
        for line in file:
            protein_id = line.strip()
            output_file = f'ipg/{protein_id}.txt'
            # Fetch IPG file using Entrez
            handle = Entrez.efetch(db="protein", id=protein_id, rettype="ipg", retmode="text")
            ipg_data = StringIO(handle.read().decode())  # Read the response data
            ipg_data_df = pd.read_csv(ipg_data, sep="\t")
            with open(output_file, 'w') as ipg_file:
                ipg_data_df.to_csv(output_file, index=False, header=False)  # Write the decoded data to the file

def create_summary_file():
    """
    Create a summary file containing the second line of each IPG file.
    """
    summary_data = []
    for filename in os.listdir('ipg/'):
        with open(os.path.join('ipg/', filename), 'r') as file:
            lines = file.readlines()
            second_line_values = lines[1].strip().split('\t')  # Assuming tab-separated values
            summary_data.append(second_line_values)
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

def extend_ipg_files_with_assembly_length():
    """
    Extend each IPG file with assembly length information.
    """
    accs = []
    for filename in os.listdir('ipg/'):
        accs.append(filename[:-4])  # Remove '.txt' extension
    
    assembly_lengths = fetch_assembly_lengths(accs)
    
    for acc in accs:
        ipg_file = f'ipg/{acc}.txt'
        if acc in assembly_lengths:
            assembly_length = assembly_lengths[acc]
            with open(ipg_file, 'a') as file:
                file.write(f'\nAssembly Length: {assembly_length}\n')

def fetch_assembly_lengths(accs):
    """
    Fetch assembly lengths for a list of assembly accessions.
    """
    assembly_lengths = {}
    for acc in accs:
        handle = Entrez.efetch(db="assembly", id=acc, rettype="docsum", retmode="xml")
        record = Entrez.read(handle)
        assembly_length = int(record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStats']['Length'])
        assembly_lengths[acc] = assembly_length
    return assembly_lengths

def download_genome_data():
    """
    Download genome data based on the list of assembly accessions.
    """
    subprocess.run(['../../ncbi/datasets', 'download', 'genome', 'accession', '--inputfile', 'assm_accs.txt', '--include', 'gff3'])

def unzip_downloaded_files():
    """
    Unzip the downloaded files.
    """
    subprocess.run(['unzip', 'ncbi_dataset.zip'])


