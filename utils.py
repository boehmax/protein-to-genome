import subprocess
import os
from Bio import Entrez

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
            ipg_data = handle.read()  # Read the response data
            with open(output_file, 'w') as ipg_file:
                ipg_file.write(ipg_data.decode())  # Write the decoded data to the file

def create_summary_file():
    """
    Create a summary file containing the second line of each IPG file.
    """
    summary_lines = []
    for filename in os.listdir('ipg/'):
        with open(os.path.join('ipg/', filename), 'r') as file:
            second_line = file.readlines()[1]
            summary_lines.append(second_line)
    with open('ipg_summary.txt', 'w') as summary_file:
        summary_file.write('\n'.join(summary_lines))

def process_summary_file():
    """
    Process the summary file to extract relevant information.
    """
    with open('ipg_summary.txt', 'r') as summary_file:
        lines = summary_file.readlines()
        accs_protein = []
        accs = []
        for line in lines:
            parts = line.split()
            acc = parts[-1]
            accs_protein.append(f'{acc} {parts[6]}')
            if acc.startswith('GC'):
                accs.append(acc)
        with open('assm_accs_protein.txt', 'w') as accs_protein_file:
            accs_protein_file.write('\n'.join(accs_protein))
        with open('assm_accs.txt', 'w') as accs_file:
            accs_file.write('\n'.join(accs))

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
    Entrez.email = "your_email@example.com"  # Enter your email here
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


