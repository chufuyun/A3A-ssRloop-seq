import pysam
import random
import os
import argparse

def generate_skewed_distribution(mean, max_value, size):
    values = []
    for _ in range(size):
        num_mutations = random.paretovariate(1)  # Generate a random number with a skewed distribution
        num_mutations = min(num_mutations, max_value)  # Limit to the max value
        values.append(int(num_mutations))
    return values

def convert_base(seq, base_from, base_to, mutations_count):
    mutated_seq = list(seq)
    positions = [i for i, base in enumerate(seq[:60]) if base == base_from]
    mutated_positions = random.sample(positions, min(mutations_count, len(positions)))
    for pos in mutated_positions:
        mutated_seq[pos] = base_to
    return ''.join(mutated_seq), mutated_positions

def convert_base_for_R2(seq, base_from, base_to, mutations_count):
    mutated_seq = list(seq)
    positions = [i for i, base in enumerate(seq[-60:]) if base == base_from]
    mutated_positions = random.sample(positions, min(mutations_count, len(positions)))
    for pos in mutated_positions:
        mutated_seq[-60 + pos] = base_to
    return ''.join(mutated_seq), mutated_positions

def convert_bam(bam_file_name, output_bam_name, output_txt_name):
    bam_file = pysam.AlignmentFile(bam_file_name, 'rb')
    output_bam = pysam.AlignmentFile(output_bam_name, 'wb', header=bam_file.header)

    with open(output_txt_name, 'w') as output_txt:
        for read in bam_file:
            # 跳过没有被映射的reads和长度小于60的·reads
            if read.is_unmapped or read.query_length < 60:
                continue  # 如果read没有被映射，跳到下一个read

            original_qual = read.qual

            # If the read is R1 (flag 99 or 83), convert C's to T's within the first 60 bp
            if read.flag in [99, 83]:
                mutations_count = generate_skewed_distribution(47, 2, 1)[0]
                read.seq, mutated_positions = convert_base(read.seq, 'C', 'T', mutations_count)
                if mutated_positions:
                    for pos in mutated_positions:
                        output_txt.write(f"{read.qname}\t{read.reference_name}\t{read.reference_start + pos+1}\n")

            # If the read is R2 (flag 147 or 163), convert C's to T's within the last 60 bp
            elif read.flag in [147, 163]:
                mutations_count = generate_skewed_distribution(47, 2, 1)[0]
                read.seq, mutated_positions = convert_base_for_R2(read.seq, 'C', 'T', mutations_count)
                if mutated_positions:
                    for pos in mutated_positions:
                        output_txt.write(f"{read.qname}\t{read.reference_name}\t{read.reference_start + read.query_length - 60 + pos+1}\n")

            read.qual = original_qual
            output_bam.write(read)

    bam_file.close()
    output_bam.close()


def extract_tags(bam_file, output_file):
    # Open the BAM file for reading
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Open the output file for writing
        with open(output_file, "w") as output_txt:
            # Iterate through each read in the BAM file
            for read in bam:
                # Skip if the read is unmapped
                if read.is_unmapped:
                    continue
                # Get the XM:Z tag value
                xm_tag = read.get_tag("XM", with_value_type=True)[0] if read.has_tag("XM") else None
                if xm_tag:
                    # Find all positions of x, m, or z
                    positions = [i for i, char in enumerate(xm_tag, start=1) if char in 'xzh']
                        # Format the required information
                    for pos in positions:
                        output_txt.write(f"{read.qname}\t{read.reference_name}\t{read.reference_start + pos}\n")
                            
                            
                            
                        
def compare_files(a_file, b_file):
    # Read the contents of the files
    with open(a_file, 'r') as f:
        a_content = set(f.readlines())

    with open(b_file, 'r') as f:
        b_content = set(f.readlines())

    # Calculate TP, FP, and FN using set operations
    TP = len(a_content & b_content)  # Intersection of a and b
    FN = len(a_content - b_content)  # Items in a but not in b
    FP = len(b_content - a_content)  # Items in b but not in a

    return TP, FP, FN

def calculate_metrics(TP, FP, FN):
    # Calculate precision, recall, and f1-score
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1_score = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    return precision, recall, f1_score

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="A script to simulate and analyze sequence data.")
    parser.add_argument('-1', dest='fq1', required=True, help='Path to the first FastQ file.')
    parser.add_argument('-2', dest='fq2', required=True, help='Path to the second FastQ file.')
    parser.add_argument('--genome', dest='genome_path',help='Path to the genome folder.', default="/home/limh/cfy/genome/hg38/")
    ## set seed
    parser.add_argument('-s', dest='seed', default=0, help='seed number')

    return parser.parse_args()

def bismark_command(input_R1, input_R2, genome_path, output_dir="."):
    """Run Bismark command."""
    command = f"bismark --bowtie2 --score_min L,0,-0.6 --un -p 5 -o {output_dir} {genome_path} -1 {input_R1} -2 {input_R2}"
    os.system(command)

def bismark_R1(input_R1,genome_path, output_dir="."):
    """Run Bismark command."""
    command = f"bismark --bowtie2 --score_min L,0,-0.6 --un -p 5 -o {output_dir} {genome_path} {input_R1} "
    os.system(command)

def bismark_R2(input_R2, genome_path, output_dir="."):
    """Run Bismark command."""
    command = f"bismark --bowtie2 --pbat --score_min L,0,-0.6 --un -p 5 -o {output_dir} {genome_path} {input_R2}"
    os.system(command)

import os

def samtools_sort(bam_file,bam_sorted_file):
    """Index a BAM file."""
    command = f"samtools sort -o {bam_sorted_file} {bam_file}"
    os.system(command)

def samtools_index(bam_file):
    """Index a BAM file."""
    command = f"samtools index {bam_file}"
    os.system(command)

def bedtools_bamtofastq(bam_file, fq1, fq2):
    """Convert BAM to FastQ."""
    command = f"bedtools bamtofastq -i {bam_file} -fq {fq1} -fq2 {fq2}"
    os.system(command)

def modify_output_file(input_file, output_file):
    """Modify output file using awk."""
    command = f'''awk -F "\t" 'BEGIN {{OFS=FS}} {{sub(/\\/[12]$/, "", $1); print}}' {input_file} > {output_file}'''
    os.system(command)

def main(args):
    # Set random seed
    s=args.seed
    random.seed(s)
    input_R1 = args.fq1
    input_R2 = args.fq2
    genome_path = args.genome_path
    # Run Bismark
    bismark_command(input_R1, input_R2, genome_path)
    
    # sort BAM
    bam_file = "test.R1_bismark_bt2_pe.bam"
    bam_sorted_file = "test.sorted.bam"
    samtools_sort(bam_file, bam_sorted_file)
    # Index BAM
    bam_file = "test.sorted.bam"
    samtools_index(bam_sorted_file)
    
    # Convert BAM and get mutation statistics
    output_bam = "test.bam"
    output_txt = "test.txt"
    convert_bam(bam_sorted_file, output_bam, output_txt)
    
    # Convert BAM to FastQ
    converted_R1 = "test.1.fq"
    converted_R2 = "test.2.fq"
    bedtools_bamtofastq(output_bam, converted_R1, converted_R2)
    
    # Run Bismark again on the converted FastQ
    bismark_command(converted_R1, converted_R2, genome_path)
    unmapped_R1 = "test.1.fq_unmapped_reads_1.fq.gz"
    bismark_R1(unmapped_R1, genome_path)
    unmapped_R2 = "test.2.fq_unmapped_reads_2.fq.gz"
    bismark_R2(unmapped_R2, genome_path)
    
    # Sort the new BAM
    new_bam_file = "test.1_bismark_bt2_pe.bam"
    new_bam_sorted_file = "test.1.sorted.bam"
    R1_bam_file = "test.1.fq_unmapped_reads_1_bismark_bt2.bam"
    R1_bam_sorted_file = "test.1.unmapped.sorted.bam"
    R2_bam_file = "test.2.fq_unmapped_reads_2_bismark_bt2.bam"
    R2_bam_sorted_file = "test.2.unmapped.sorted.bam"


    samtools_sort(new_bam_file, new_bam_sorted_file)
    samtools_sort(R1_bam_file,R1_bam_sorted_file)
    samtools_sort(R2_bam_file,R2_bam_sorted_file)
    # Index the new BAM
    samtools_index(new_bam_sorted_file)
    samtools_index(R1_bam_sorted_file)
    samtools_index(R2_bam_sorted_file)

    # Extract tags
    tags_output = "PE.txt"
    extract_tags(new_bam_sorted_file, tags_output)
    # Modify output file
    modified_output = "PE_output.txt"
    modify_output_file(tags_output, modified_output)
    
    # Compare files and calculate metrics for PE
    TP, FP, FN = compare_files(output_txt, modified_output)
    precision, recall, f1_score = calculate_metrics(TP, FP, FN)
    # Output results
    with open('precision_PE.txt', 'w') as f:
        f.write(f"Precision: {precision:.4f}\n")
        f.write(f"Recall: {recall:.4f}\n")
        f.write(f"F1-Score: {f1_score:.4f}\n")
    # Compare files and calculate metrics for SE
    tags_output = "SE.R1.txt"
    extract_tags(R1_bam_sorted_file, tags_output)
    modified_output = "SE.R1_output.txt"
    modify_output_file(tags_output, modified_output)
    tags_output = "SE.R2.txt"
    extract_tags(R2_bam_sorted_file, tags_output)
    modified_output = "SE.R2_output.txt"
    modify_output_file(tags_output, modified_output)
    os.system('cat SE.R1.txt SE.R2_output.txt PE_output.txt > SE_output.txt ')
    modified_output = "SE_output.txt"
    TP, FP, FN = compare_files(output_txt, modified_output)
    precision, recall, f1_score = calculate_metrics(TP, FP, FN)
    # Output results
    with open('precision_SE.txt', 'w') as f:
        f.write(f"Precision: {precision:.4f}\n")
        f.write(f"Recall: {recall:.4f}\n")
        f.write(f"F1-Score: {f1_score:.4f}\n")

if __name__ == "__main__":
    args = parse_args()
    main(args)
