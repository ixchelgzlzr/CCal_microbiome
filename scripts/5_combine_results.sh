# combine metaphlan output

# go to dir with all the profiles
cd output/metaphlan_out/profiles

# Remove metaphlan from the name of a file 
for f in *_metaphlan.txt; do
    mv "$f" "${f/_metaphlan/}"
done


# then merge all of them in a single file
merge_metaphlan_tables.py *.txt > merged_abundance_table.txt

# trying some of their scripts
mkdir tests
sgb_to_gtdb_profile.py -i merged_abundance_table.txt -o tests/metaphlan_output_gtdb.txt