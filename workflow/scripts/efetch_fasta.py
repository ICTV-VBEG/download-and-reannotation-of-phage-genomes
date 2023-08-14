#!/usr/bin/env python3
from Bio import Entrez

def download_fasta_from_ncbi(accession_file, output_file):
    """
    Downloads FASTA files from NCBI given a file with a list of accession IDs.

    Parameters:
        accession_file (str): Path to the file containing a list of accession IDs, one per line.
        output_directory (str): Directory where the downloaded FASTA files will be saved.

    Returns:
        None
    """
    # Read the accession IDs from the file
    print(accession_file,output_file)
    with open(accession_file, 'r') as file:
        accession_ids = [line.strip() for line in file]

    # Use the Entrez module to fetch FASTA sequences from NCBI
    Entrez.email = 'example@example.com'  # this should come from the configfile

    for accession_id in accession_ids:
        try:
            handle = Entrez.efetch(db='nucleotide', id=accession_id, rettype='fasta', retmode='text')
            fasta_data = handle.read()
            handle.close()

            with open(output_file, 'a') as fasta_file:
                fasta_file.write(fasta_data)

            print(f'Successfully downloaded {accession_id}.')
        except Exception as e:
            print(f'Error downloading {accession_id}: {str(e)}')


# Example usage:
download_fasta_from_ncbi(snakemake.input[0], snakemake.output[0])
