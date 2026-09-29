#!/usr/bin/env python3
"""
统一所有物种的FASTA header格式
格式: >物种名|蛋白质ID
"""
import sys
import gzip
import re


def read_fasta(file_path):
    """读取FASTA文件，返回{header: sequence}"""
    sequences = {}
    current_header = None
    
    open_func = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'
    
    with open_func(file_path, mode) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_header = line[1:]
                sequences[current_header] = []
            elif current_header:
                sequences[current_header].append(line)
    
    return {k: ''.join(v) for k, v in sequences.items()}


def extract_protein_id(header, species_name):
    """从header中提取蛋白质ID"""
    # 格式: >物种名|蛋白质ID 或 >蛋白质ID 或 >物种名|其他|蛋白质ID
    
    # 情况1: 已经是我们想要的格式 >物种名|蛋白质ID
    if f">{species_name}|" in f">{header}|":
        parts = header.split('|')
        if len(parts) >= 2:
            return parts[-1]  # 取最后一部分作为蛋白质ID
    
    # 情况2: >Acch|PSS32443.1 -> PSS32443.1
    if '|' in header:
        return header.split('|')[-1]
    
    # 情况3: 直接是蛋白质ID
    return header.split()[0]


def standardize_header(input_file, output_file, species_name):
    """标准化header格式"""
    print(f"读取: {input_file}")
    sequences = read_fasta(input_file)
    print(f"找到 {len(sequences)} 个序列")
    
    # 写入标准化后的文件
    with open(output_file, 'w') as f:
        count = 0
        for header in sorted(sequences.keys()):
            protein_id = extract_protein_id(header, species_name)
            if protein_id:
                f.write(f">{species_name}|{protein_id}\n")
                seq = sequences[header]
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i + 80] + '\n')
                count += 1
    
    print(f"输出: {output_file} ({count} 个序列)")
    return count


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python standardize_headers.py <input.fa> <output.fa> <species_name>")
        sys.exit(1)
    
    standardize_header(sys.argv[1], sys.argv[2], sys.argv[3])
