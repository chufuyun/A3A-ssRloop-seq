# -*- coding: utf-8 -*-
import pysam
import argparse
import time
import os

# 解析命令行参数
parser = argparse.ArgumentParser(description='Get the bam file from the command line')
parser.add_argument('-i', '--input', help='input bam file', required=True)
parser.add_argument('-o', '--output', help='output bam file id', required=True)
parser.add_argument('-t1', '--threshold1', help='the minimum number of C>T sites', default=1)
parser.add_argument('-t2', '--threshold2', help='the minimum number of C content')
parser.add_argument('-t3', '--threshold3', help='the maximum number of C content')
parser.add_argument('-s', '--summary', help='the summary file', default='summary.txt')

args = parser.parse_args()

# 获取开始时间
start_time = time.time()
#records = []

class_folders = ['ClassI', 'ClassII', 'Unclassified']
for folder in class_folders:
    try:
        os.makedirs(folder)
    except FileExistsError:
        pass

t1 = args.threshold1
t2 = args.threshold2
t3 = args.threshold3

# 定义分类函数
def classify_read(t2, t3):
    if float(t2) >= 40 and float(t3) <= 200:
        return 'ClassI'
    elif float(t2) >= 0 and float(t3) <= 40:
        return 'ClassII'
    else:
        return 'Unclassified'


# 定义处理read的函数
def process_read(read, class_bamfiles, t2, t3):
    xm_tag_value = read.get_tag("XM")

    if read.flag == 99 or read.flag == 83:
        if xm_tag_value is not None:
            unmethylated = xm_tag_value.count("x") + xm_tag_value.count("h") + xm_tag_value.count("z")
            methylated = xm_tag_value.count("X") + xm_tag_value.count("H") + xm_tag_value.count("Z")
            c_count = unmethylated + methylated
            g_count = read.query_alignment_sequence.count("G")
    elif read.flag == 147 or read.flag == 163:
        if xm_tag_value is not None:
            unmethylated = xm_tag_value.count("x") + xm_tag_value.count("h") + xm_tag_value.count("z")
            methylated = xm_tag_value.count("X") + xm_tag_value.count("H") + xm_tag_value.count("Z")
            c_count = unmethylated + methylated
            g_count = read.query_alignment_sequence.count("C")

    class_name = classify_read(t2, t3)
    output_bamfile_name = None  # Initialize to None

    if unmethylated >= float(t1):
        #records.append((unmethylated, methylated))
        if (read.flag == 99 or read.flag == 147) and c_count >= float(t2) and c_count <= float(t3):
            class_folder = class_name
            output_bamfile_name = os.path.join(
                class_folder,
                f"{args.output}_{class_name}_fwd_flt.bam"
            )
        elif (read.flag == 83 or read.flag == 163) and c_count >= float(t2) and c_count <= float(t3):
            class_folder = class_name
            output_bamfile_name = os.path.join(
                class_folder,
                f"{args.output}_{class_name}_rev_flt.bam"
            )
    if output_bamfile_name:
        if output_bamfile_name not in class_bamfiles[class_folder]:
            class_bamfiles[class_folder][output_bamfile_name] = pysam.AlignmentFile(output_bamfile_name, "wb", template=bamfile)

        class_bamfiles[class_folder][output_bamfile_name].write(read)

# 打开输入BAM文件
bamfile = pysam.AlignmentFile(args.input, "rb")

class_bamfiles = {folder: {} for folder in class_folders}

# 检查是否未提供阈值，如果是，则使用默认阈值列表 (t2=40, t3=200) 还有（t2=0,t3=40）
if t2 is None or t3 is None:
    for read in bamfile:
        process_read(read, class_bamfiles, 40, 200)
else:
    for read in bamfile:
        process_read(read, class_bamfiles, t2, t3)

bamfile.close()

# 初始化
class_bamfiles = {folder: {} for folder in class_folders}

for class_folder, bamfiles in class_bamfiles.items():
    for output_bamfile in bamfiles.values():
        output_bamfile.close()

# 打印总读数和程序运行时间
end_time = time.time()
print("Total run time:", end_time - start_time)

# with open(args.summary, "w") as f:
#     f.write("Unmethylated_C\tMethylated_C\n")  # 列名
#     for record in records:
#         f.write(f"{record[0]}\t{record[1]}\n")
