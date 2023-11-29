## write a bash script for bismark alignment
## Step 1: PE alignment

#!/bin/bash

# genome_dir
genome_dir_mm10="/home/limh/cfy/genome/mm10/bwa_index/"
genome_dir_hg38="/home/limh/cfy/genome/hg38/"

if [ $# -ne 3 ]; then
    echo "usage: $0 <mm10 or hg38> <input_dir> <output_dir>"
    exit 1
fi

if [ "$1" == "mm10" ]; then
    genome_dir="$genome_dir_mm10"
elif [ "$1" == "hg38" ]; then
    genome_dir="$genome_dir_hg38"
else
    echo "请提供有效的基因组参数（mm10或hg38）"
    exit 1
fi

input_dir="$2" # reads are look like ${id}.R2.trimmed.fastq.gz
output_dir="$3"

for file in "$input_dir"/*.R1.trimmed.fastq.gz; do
    sample_name=$(basename "$file" | cut -d '.' -f 1)
    R1_file="$file"
    R2_file="${file/.R1.trimmed.fastq.gz/.R2.trimmed.fastq.gz}"


bismark --bowtie2 -N 0 -p 2 -o "$output_dir/" $genome_dir -1 "$R1_file" -2 "$R2_file"
mv "$output_dir/$sample_name*pe.bam" "$output_dir/$sample.pe.bam"
deduplicate_bismark -p --output_dir "$output_dir/" "$output_dir/$sample_name.pe.bam"
samtools sort -@ 5 -o "$output_dir/$sample.sorted.bam""$output_dir/$sample.pe.deduplicated.bam"


done

