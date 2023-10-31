#!/bin/bash

# Function to call peaks
callpeak() {
    local treat_files=()
    local control_file=""
    local s=""

    while [[ $# -gt 0 ]]; do
        case "$1" in
            -t)
                shift
                while [[ $# -gt 0 && ! $1 == -* ]]; do
                    treat_files+=("$1")
                    shift
                done
                ;;
            -c)
                shift
                control_file="$1"
                shift
                ;;
            -s)
                shift
                s="$1"
                shift
                ;;
        esac
    done

    # Ensure at least one treat file and control file is provided
    if [ -z "$control_file" ] || [ ${#treat_files[@]} -eq 0 ]; then
        echo "Usage: ./your_script.sh -t <treat1.bam> [<treat2.bam> ...] -c <control.bam> -s <species>"
        exit 1
    fi

    source /home/limh/miniconda3/bin/activate KAS-pipe2

    if [ "$s" == "mm10" ]; then
        bed="/home/limh/cfy/KAS-pipe2/annotation/mm10/mm10_Refseq.bed"
        blacklist="/home/limh/cfy/KAS-pipe2/blacklist/mm10-blacklist.bed"
        size="/home/limh/cfy/KAS-pipe2/chrom_size/mm10.chrom.sizes"
        genome_size="2707406244"
        g="mm"
    elif [ "$s" == "hg38" ]; then
        bed="/home/limh/cfy/KAS-pipe2/annotation/hg38/hg38_Refseq.bed"
        blacklist="/home/limh/cfy/KAS-pipe2/blacklist/hg38-blacklist.bed"
        size="/home/limh/cfy/KAS-pipe2/chrom_size/hg38.chrom.sizes"
        genome_size="3113504866"
        g="hs"
    fi

    # Calculate bamCoverage for control BAM file
    all_files=("$control_file" "${treat_files[@]}")
    control_id=$(basename $control_file .bam)
    for bam_file in "${all_files[@]}"; do
        id=$(basename $bam_file .bam)
        samtools index $bam_file

        # Calculate bamCoverage for BAM files
        bamCoverage -b $bam_file --outFileFormat bedgraph -p 10 -bl $blacklist --effectiveGenomeSize $genome_size --normalizeUsing RPKM -o ${id}.RPKM.bg
        macs2 callpeak -B -t ${id}.RPKM.bg -c ${control_id}.RPKM.bg -n ${id} -g ${g} -q 0.05 --format BED --nomodel

        macs2 bdgcmp -t ${id}_treat_pileup.bdg -c ${id}_control_lambda.bdg -o ${id}_FE.bdg -m FE
        sort -k1,1 -k2,2n ${id}_FE.bdg -o ${id}.FE.sort.bdg
        bedClip ${id}.FE.sort.bdg ${size} ${id}.FE.sort.output
        bedGraphToBigWig ${id}.FE.sort.output ${size} ${id}.bw
    done

    # Create an array of .bw files in the current directory
    bw_files=()
    for file in *.bw; do
        [ -e "$file" ] || continue  # Check if the file exists
        bw_files+=("$file")
    done

    computeMatrix scale-regions -S "${bw_files[@]}" -R $bed -p 15 -b 3000 -a 3000 --regionBodyLength 6000 --skipZeros --missingDataAsZero -o Rloop.mat.gz
    plotProfile -m Rloop.mat.gz -out Rloop.pdf --dpi 300 --legendLocation upper-right -y Frequency --plotHeight 7 --plotWidth 9 --perGroup
}

# Call the function with command line arguments
callpeak "$@"





