import argparse
import pysam
from collections import defaultdict

def main(args):
    # 打开输入的BAM文件
    input_bam = pysam.AlignmentFile(args.input, "rb")

    # 创建输出的BAM文件
    output_bam = pysam.AlignmentFile(args.output, "wb", template=input_bam)

    coverage_dict = defaultdict(lambda: defaultdict(int))

    # 遍历BAM文件的每个read
    for read in input_bam:
        tags = dict(read.tags)
        if "XM" in tags:
            xm = tags["XM"]
            position = read.reference_start
            strand = '+' if read.flag in {99, 83} else '-'  # 确定链方向
            for i, char in enumerate(xm):
                if char != '.':  # Ignore non-C bases
                    site_position_strand = f"{read.reference_name}:{position + i}:{strand}"  # 创建包含位点和链方向的键
                    if char.islower():  # 检查字符是否为小写，表示转换
                        coverage_dict[site_position_strand]['lower'] += 1
                    else:
                        coverage_dict[site_position_strand]['upper'] += 1

    # 输出位点的统计信息到文本文件
    with open(args.summary, "w") as f:
        f.write("chr\tpos\tstrand\tuppercase\tlowercase\tconversion rate\n")
        for site_position_strand, counts in coverage_dict.items():
            chr, pos, strand = site_position_strand.split(':')
            uppercase = counts['upper']
            lowercase = counts['lower']
            total = uppercase + lowercase
            conversion_rate = lowercase / total if total > 0 else 0
            f.write(f"{chr}\t{pos}\t{strand}\t{uppercase}\t{lowercase}\t{conversion_rate:.2f}\n")

    # 根据覆盖率和转换率阈值筛选reads
    # input_bam.reset()  # 重置BAM文件的读取指针

    # for read in input_bam:
    #     tags = dict(read.tags)
    #     if "XM" in tags:
    #         xm = tags["XM"]
    #         position = read.reference_start
    #         strand = '+' if read.flag in {99, 83} else '-'  # 确定链方向
    #         for i, char in enumerate(xm):
    #             site_position_strand = f"{read.reference_name}:{position + i}:{strand}"  # 创建包含位点和链方向的键
    #             if char.islower():
    #                 coverage = coverage_dict[site_position_strand]['lower'] + coverage_dict[site_position_strand]['upper']
    #                 total_conversions = coverage_dict[site_position_strand]['lower']
    #                 conversion_rate = total_conversions / coverage if coverage > 0 else 0
    #                 if coverage > 1 and conversion_rate > 0.01:
    #                     output_bam.write(read)
    #                     break  # 如果read满足条件，写入到输出BAM文件并跳出循环

    # input_bam.close()
    # output_bam.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process BAM file to generate statistics and filtered output.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input BAM file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output BAM file.')
    parser.add_argument('-s', '--summary', default='summary.txt', help='Path to the output summary text file.')
    args = parser.parse_args()
    main(args)
