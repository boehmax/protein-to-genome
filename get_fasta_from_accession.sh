#!/bin/bash

# Check if the input file is provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 input_file output_file"
  exit 1
fi

input_file="$1"
output_file="$2"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "Error: File '$input_file' not found."
  exit 1
fi

# Clear the output file if it exists
> "$output_file"

# Read each line from the input file and feed it to the command
while IFS= read -r line; do
  # Replace 'your_command' with the actual command you want to run
  # For example, if you want to use 'efetch' with each line as an ID:
  efetch -db nuccore -id "$line" -format fasta >> "$output_file"
done < "$input_file"