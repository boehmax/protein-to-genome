from utils import * 
from Bio import Entrez
#PIGI = Protein I gave as Input (accession number)

def main():
    print("Welcome to the Genome Data Downloader!")
    print("This program will download genome data for a list of protein accessions.")
    print("Please make sure you have the 'datasets' and 'dataformat' command line tool installed from NCBI Datasets.")
    print("You can find the installation instructions here: https://www.ncbi.nlm.nih.gov/datasets/docs/command-line/")
    print("Please also make sure you have an API key from NCBI. You can get one here: https://www.ncbi.nlm.nih.gov/account/settings/")
    print("\n")
    print("Let's get started!")
    # Define the input file containing protein IDs
    filename = input("Give path to your file of protein accessions (default: input/example_proteins.txt): ")
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
    
    print("Genome data download complete!")
    print("You can find the downloaded files in the 'ncbi/datasets' directory.")
    print("Thank you for using the Genome Data Downloader! :)")

if __name__ == "__main__":
    main()
