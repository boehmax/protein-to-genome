# Protein to Genome 

This project focuses on connecting proteins to their corresponding genome assembly in a batch fasion. It is designed to assist bioinformatics researchers in understanding the relationship between proteins and their genetic origins.

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Additional Requirements](#additional_requirements)
## Introduction
The `protein-to-genome` project provides tools and scripts to map protein accesions to genom assemblies. This can be useful for various bioinformatics analyses, including gene prediction, annotation, and evolutionary studies.

## Features
- Collects all Identical Protein accessions
- Chooses Genome Assembly with best stats
- Downloads all genomes
- Exports summary files for helpfull for further prosessing (e.g. boehmax/protein-neighbours)

## Installation
To install the necessary dependencies, run:
```bash
pip install -r requirements.txt
```
## Additional Requirements

Please make sure you have the `datasets` and `dataformat` command-line tools installed from NCBI Datasets. You can find the installation instructions here:
https://www.ncbi.nlm.nih.gov/datasets/docs/command-line/ . For smooth sailing best to install it into your homedirectory.

Please also make sure you have an API key from NCBI. You can get one here:
https://www.ncbi.nlm.nih.gov/account/settings/

## Usage
To map a list of protein accession to their genomes, use the following command:
```bash
python3 main.py -f <protein_acccesions_file> -e <Your.email@from.ncbi.account.com> -k <your_NCBI_API_key> 
```

### Output

#### IPG folder
After running the script, you will get a folder containing `.csv` files. Each file is named after the protein accession provided as input and is stored in the following directory structure: output\current_date\ipg
Where `current_date` is the date when the script was run.

Each CSV file contains the following columns:
- `accver`: Accession version of the protein you gave as input.
- `name`: Name of the protein.
- `taxid`: Taxonomic ID of the organism.
- `slen`: Sequence length of the protein.
- `org`: Organism name.
- `kingdom_taxid`: Taxonomic ID of the kingdom.
- `kingdom`: Kingdom name.
- `accver_dup`: Accession number of the identical proteins found.
- `source`: Source of the data.
- `name_dup`: Name of the protein of the identical proteins found.
- `taxid_dup`: Duplicate taxonomic ID.
- `org_dup`: Duplicate organism name.
- `kingdom_taxid_dup`: Duplicate kingdom taxonomic ID.
- `kingdom_dup`: Duplicate kingdom name.
- `priority`: Priority of the data.
- `accver_dup_dup`: 
- `start`: Start position of the sequence.
- `stop`: Stop position of the sequence.
- `strand`: Strand information.
- `taxid_dup_dup`: Second duplicate taxonomic ID.
- `org_dup_dup`: Second duplicate organism name.
- `kingdom_taxid_dup_dup`: Second duplicate kingdom taxonomic ID.
- `kingdom_dup_dup`: Second duplicate kingdom name.
- `strain`: Strain information.
- `assembly`: Assembly information.
- `PIGI`: Protein i gave as input, should be the same as `accver`.
- `assembly_accession`: Assembly accession number.
- `checkm_completeness`: Completeness score from CheckM.
- `checkm_contamination`: Contamination score from CheckM.
- `checkm_version`: Version of CheckM used.
- `contig_N50`: N50 value of contigs.
- `contig_L50`: L50 value of contigs.
- `ungaped_seq_len`: Length of the ungapped sequence.
- `seq_len`: Total sequence length.

#### Summary Folder
In addition to the `.csv` files, you will also get a summary folder (output\current_date\summary) containing the following files:
- `assm_accs.csv`: A list of all genome assembly numbers downloaded.
- `assm_accs_protein.csv`: A table containing all genomes downloaded together with the protein you are interested in.
- `ipg_summary.csv`: A summary file with all proteins and their genome assembly in the same structure as the `.csv` files in the `ipg` folder. These were the ones that were picked based on the best genome assembly.
- `ncbi_dataset`: The NCBI dataset that was downloaded and unzipped, containing all `.gff` genome files.
- `ncbi_dataset.zip`: The zipped version of the `ncbi_dataset` folder in case of limited storage.
- `README.md`: A README file provided by NCBI to explain the dataset that was just downloaded.

This detailed output provides comprehensive information about each protein accession and its corresponding genomic data, which can be useful for further bioinformatics analyses.

## Contributing
Contributions are welcome! Please fork the repository and submit a pull request. Understand that I am primarily a wet lab biochemist and it might take some time, until I can consider your contribution. Thank you!

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE.md) file for details.