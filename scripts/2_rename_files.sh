#!/usr/bin/env bash

match_data=data/matching_table_batch2.csv
in_dir=data/clean_reads
out_dir=data/renamed_reads

# make outout directory
mkdir -p "$out_dir"

# In this script, I iterate over the rows on the file samples_names.txt
# locate the file 
while IFS=, read -r new_basename r1_file r2_file; do

    # build file names
    this_R1_in="${in_dir}/${r1_file%$'\r'}"
    this_R2_in="${in_dir}/${r2_file%$'\r'}"

    this_R1_out="${out_dir}/${new_basename%$'\r'}_R1.fastq.gz"
    this_R2_out="${out_dir}/${new_basename%$'\r'}_R2.fastq.gz"

    # If file exist
    if [[ -f "$this_R1_in" && -f "$this_R2_in" ]]; then

    # copy and rename in a new folder
    cp -v "$this_R1_in" "$this_R1_out"
    cp -v "$this_R2_in" "$this_R2_out"

    else

    # throw a warning
    echo "WARNING: missing input for mapping: new='$new_basename'" >&2
    [[ -f "$this_R1_in" ]] || echo "  missing: $this_R1_in" >&2
    [[ -f "$this_R2_in" ]] || echo "  missing: $this_R2_in" >&2

    fi

done < "$match_data"