#!/usr/bin/env python3
import sys

def split_multi_fasta(input_file, output_paths):
    """
    Splits a multi FASTA file into separate files and saves each sequence into specified output paths.

    Parameters:
        input_file (str): Path to the input multi FASTA file.
        output_paths (list): List of output file paths where individual sequences will be saved.

    Returns:
        None
    """

    def write_fasta_record(output_path, header, sequence):
        with open(output_path, 'a') as output_file:
            output_file.write(f'{header}\n{sequence}\n')

    current_output_path_index = 0
    current_output_path = output_paths[current_output_path_index]

    with open(input_file, 'r') as input_file:
        header = ''
        sequence = ''

        for line in input_file:
            line = line.strip()

            if line.startswith('>'):
                if header:
                    write_fasta_record(current_output_path, header, sequence)
                    sequence = ''

                header = line
                current_output_path = output_paths[current_output_path_index]
                current_output_path_index = (current_output_path_index + 1) % len(output_paths)

            else:
                sequence += line

        # Write the last sequence after the loop ends
        if header:
            write_fasta_record(current_output_path, header, sequence)


if __name__ == "__main__":
    # Check if the correct number of arguments are provided
    if len(sys.argv) < 3:
        print("Usage: python script_name.py input_file output_path1 output_path2 ...")
        sys.exit(1)

    # First argument is the input file path
    input_fasta_file = sys.argv[1]

    # Remaining arguments are the output file paths
    output_paths_list = sys.argv[2:]

    # Call the function to split the multi FASTA file
    split_multi_fasta(input_fasta_file, output_paths_list)

