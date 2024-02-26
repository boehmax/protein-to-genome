from utils import * 
from Bio import Entrez


def main():
    # Define the input file containing protein IDs
    filename = input("Give path to your file of protein accessions: ")
    if filename == '':
        # If it is, use a default value
        filename = 'input/example_proteins.txt'
    Entrez.email = input("Your Email Adress: ") 
    Entrez.api_key = input("Your API Key: ")
    # Step 1: Retrieve IPG files for each protein ID
    retrieve_protein_info(filename)
    
    # Step 2: Create a summary file containing the second line of each IPG file
    create_summary_file()
    
    # Step 3: Extend IPG files with assembly length information WORK IN PROGESS!!
    # extend_ipg_files_with_assembly_length()
    
    # Step 4: Process the summary file to extract relevant information
    process_summary_file()
    
    # Step 5: Download genome data based on the list of assembly accessions
    download_genome_data()
    
    # Step 6: Unzip the downloaded files
    unzip_downloaded_files()

if __name__ == "__main__":
    main()
