from utils import * 
from Bio import Entrez
import argparse
#PIGI = Protein I gave as Input (accession number)

def main():
    """
    Main function to download genome data for a list of protein accessions.
    It retrieves IPG files, extends them with assembly information, and generates a summary file.
    After processing the summary file, it downloads genome data for the best genome assembly for each protein.
    """
    print("Welcome to the Genome Data Downloader!")
    print("This program will download genome data for a list of protein accessions.")
    print("Please make sure you have the 'datasets' and 'dataformat' command line tool installed from NCBI Datasets.")
    print("You can find the installation instructions here: https://www.ncbi.nlm.nih.gov/datasets/docs/command-line/")
    print("Please also make sure you have an API key from NCBI. You can get one here: https://www.ncbi.nlm.nih.gov/account/settings/")
    print("\n")
    print("Let's get started!")
    # Define the input file containing protein IDs
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Download genome data for protein accessions.")
    parser.add_argument('-f', '--file', type=str, default='input/example_proteins.txt', help='Path to your file of protein accessions (default: input/example_proteins.txt)')
    parser.add_argument('-e', '--email', type=str, required=True, help='Your Email Address')
    parser.add_argument('-k', '--api_key', type=str, required=True, help='Your NCBI API Key')
    args = parser.parse_args()
    # Use the parsed arguments
    filename = args.file
    Entrez.email = args.email
    Entrez.api_key = args.api_key
    
    # Step 1.1: Retrieve IPG files for each protein ID
    retrieve_protein_info(filename)
    
    # Step 1.2: Extend IPG files with assembly information and sorts them, so that best assembly is first
    extend_ipg_files_with_assembly_information(api_key=args.api_key)
    
    # Step 1.3: Generate a file containing the protein alias
    generate_protein_alias_ipg()
    
    # Step 2: Create a summary file containing the second line of each IPG file
    create_summary_file()
        
    # Step 3: Process the summary file to extract relevant information
    process_summary_file()
    
    # Step 4: Download genome data based on the list of assembly accessions
    download_genome_data(api_key=args.api_key)
    
    # Step 5: Unzip the downloaded files
    unzip_downloaded_files()
    
    print("Genome data download complete!")
    print(f"You can find the downloaded files in the 'output/YYYY-MM-DD/summary/ncbi/datasets' directory.")
    print("Thank you for using the Genome Data Downloader! :)")

if __name__ == "__main__":
    main()
