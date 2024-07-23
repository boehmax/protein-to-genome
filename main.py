from utils import * 
from Bio import Entrez
#PIGI = Protein I gave as Input (accession number)

def main():
    # Define the input file containing protein IDs
    filename = input("Give path to your file of protein accessions: ")
    if filename == '':
        # If it is, use a default value
        filename = 'input/example_proteins.txt'
    Entrez.email = input("Your Email Adress: ") 
    Entrez.api_key = input("Your API Key: ")
    # Step 1.1: Retrieve IPG files for each protein ID
    retrieve_protein_info(filename)
    
    # Step 1.2: Extend IPG files with assembly information and sorts them, so that best assembly is first
    extend_ipg_files_with_assembly_information()
    
    # Step 1.3: Generate a file containing the protein alias
    generate_protein_alias_ipg()
    
    # Step 2: Create a summary file containing the second line of each IPG file
    create_summary_file()
        
    # Step 3: Process the summary file to extract relevant information
    process_summary_file()
    
    # Step 4: Download genome data based on the list of assembly accessions
    download_genome_data()
    
    # Step 5: Unzip the downloaded files
    unzip_downloaded_files()

if __name__ == "__main__":
    main()
