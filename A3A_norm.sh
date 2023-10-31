#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <bam_file1> <bam_file2> [<bam_file3> ...] <annotation> <fwd/rev>"
    echo "Annotation options: hg38, mm10"
    exit 1
fi

# Get the annotation and parameter arguments
annotation="${@: -2:1}"  # Get the second to last argument
param="${@: -1}"       # Get the last argument

# Output the annotation and parameter for debugging
echo "Annotation: $annotation"
echo "Parameter: $param"

# Check if annotation is valid
if [ "$annotation" != "hg38" ] && [ "$annotation" != "mm10" ]; then
    echo "Invalid annotation '$annotation'. Use 'hg38' or 'mm10'."
    exit 1
fi

# Determine the BED file, blacklist file, and effectiveGenomeSize based on the annotation
if [ "$annotation" == "hg38" ]; then
    bed="/home/limh/cfy/KAS-pipe2/annotation/hg38/hg38_Refseq.bed"
    blacklist="/home/limh/cfy/KAS-pipe2/blacklist/hg38-blacklist.bed"
    genome_size="3113504866"
elif [ "$annotation" == "mm10" ]; then
    bed="/home/limh/cfy/KAS-pipe2/annotation/mm10/mm10_Refseq_${param}.bed"
    blacklist="/home/limh/cfy/KAS-pipe2/blacklist/mm10-blacklist.bed"
    genome_size="2707406244"
fi

# Loop through all input BAM files
for ((i = 1; i <= $# - 2; i++)); do
    bam_file="${!i}"
    id=$(basename "$bam_file" .bam)

    # Index the BAM file
    samtools index "$bam_file"

    source /home/limh/miniconda3/bin/activate KAS-pipe2

    # Normalize the BAM file
    bamCoverage -b "$bam_file" -o "${id}.bw" -p 10 -bl "$blacklist" --effectiveGenomeSize "$genome_size" --normalizeUsing RPKM
done

# Find all ${param}.bw files and add them to bw_files array
bw_files=()
for file in *${param}.bw; do
    [ -e "$file" ] || continue  # Check if the file exists
    bw_files+=("$file")
done

# Compute the matrix for all BAM files
computeMatrix scale-regions -S "${bw_files[@]}" -R "$bed" -p 15 -b 3000 -a 3000 --regionBodyLength 6000 --skipZeros --missingDataAsZero -o "${param}.mat.gz"

# Plot the profile
plotProfile -m "${param}.mat.gz" -out "${param}.pdf" --dpi 300 --legendLocation upper-right -y Frequency --plotHeight 7 --plotWidth 9 --perGroup
