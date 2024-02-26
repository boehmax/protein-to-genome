# From Protein to Genome

This script performs a series of operations on protein data. It retrieves IPG files for each protein ID listed in an input file, creates a summary file, processes the summary file to extract relevant information, downloads genome data based on the list of assembly accessions, and unzips the downloaded files.

## Prerequisites

The script requires the following Python packages:

- Biopython
- pandas
- numpy
- os
- subprocess

You can install these packages using pip:

```bash
pip install biopython pandas numpy
```

The script requires the following additional software:

- [NCBI command lines tools](https://github.com/ncbi/datasets): Please download them and add the path to the function download_genome_data (utils.py) in the PATH_TO_NCBI_DATASETS variable.
- *unzip*

You can install *unzip* using:

```bash
sudo apt install unzip
```

It is also recommended to use the NCBI databases with an account that you can create at their webpage. An **API key** is not necessary to enter for this script to run but if you are a heavy user, please get an API key and use it in your request. Lern how, [here](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/).

## Usage
To run the script, use the following command:

```bash
python main.py
```

When prompted, enter the path to your file of protein accessions. If you don't enter anything, the script will use 'input/example_proteins.txt' as a default.

Next, enter your email address and API key for NCBI Entrez.

The script will then perform the following steps:

1. Retrieve IPG files for each protein ID listed in the input file.
2. Create a summary file containing the second line of each IPG file.
3. (Work in progress) Extend IPG files with assembly length information. And sorts them by lenght in a descending order, to get the most compelte assembly for a protein.
4. Process the summary file to extract genome assembly accessions.
5. Download genome data based on the list of assembly accessions.
6. Unzip the downloaded files.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
