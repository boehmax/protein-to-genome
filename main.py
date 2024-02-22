from utils import * 
import subprocess
import os
from Bio import Entrez


def main():
    # Define the input file containing protein IDs
    filename = "input/example_proteins.txt"
    Entrez.email = 'maximilian.bohm@kemi.uu.se' 
    # Step 1: Retrieve IPG files for each protein ID
    retrieve_protein_info(filename)
    
    # Step 2: Create a summary file containing the second line of each IPG file
    #create_summary_file()
    
    # Step 3: Process the summary file to extract relevant information
   # process_summary_file()
    
    # Step 4: Extend IPG files with assembly length information
   # extend_ipg_files_with_assembly_length()
    
    # Step 5: Download genome data based on the list of assembly accessions
    #download_genome_data()
    
    # Step 6: Unzip the downloaded files
    #unzip_downloaded_files()

if __name__ == "__main__":
    main()
