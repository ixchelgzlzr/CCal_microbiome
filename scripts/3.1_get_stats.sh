# GET THE STATS BEFORE AND AFTER FILTERING READS

# Before filtering
seqkit stats data/renamed_reads/*.fastq.gz > before_filtering_read.stats.txt

# move unmatched reads to new directory
cd output/micro_fastq
mkdir single_reads
mv *_single.fastq.gz single_reads/
cd ..
cd ..

# After filtering
seqkit stats output/micro_fastq/*.fastq.gz > after_filtering_read.stats.txt