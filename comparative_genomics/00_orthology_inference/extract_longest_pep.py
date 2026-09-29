#!/usr/bin/env python3
"""
从蛋白质FASTA文件中提取每个基因的最长转录本
处理格式: >species|gene_id|transcript_id 或 >gene_id|transcript_id
"""
import sys
import re
from collections import defaultdict

def extract_gene_id(header):
    """从header中提取基因ID（去掉最后的转录本版本号）"""
    # 格式1: >18Diospyros_kaki|Dika|evm.model.contig20580.2 -> evm.model.contig20580
    # 格式2: >Dika|evm.model.contig20580.2 -> evm.model.contig20580
    # 格式3: >evm.model.contig20580.2 -> evm.model.contig20580
    
    parts = header.split('|')
    if len(parts) >= 3:
        # 取第三部分，去掉版本号
        gene = parts[2]
    elif len(parts) == 2:
        gene = parts[1]
    else:
        gene = header.split()[0]
    
    # 去掉最后的版本号（如.2）
    gene_base = gene.rsplit('.', 1)
    if len(gene_base) == 2 and gene_base[1].isdigit():
        return gene_base[0]
    return gene

def get_sequence_length(seq_lines):
    """计算蛋白质序列长度（去除*）"""
    seq = ''.join(seq_lines).replace('*', '').strip()
    return len(seq)

def extract_longest_transcripts(fasta_file, output_file, species_name=None):
    """从FASTA文件中提取每个基因的最长转录本"""
    gene_sequences = defaultdict(list)  # gene_id -> list of (header, seq_lines, length)
    current_gene = None
    current_seq = []
    current_header = None
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # 保存之前的序列
                if current_gene is not None:
                    gene_sequences[current_gene].append({
                        'header': current_header,
                        'seq': current_seq
                    })
                
                # 解析新序列
                current_header = line[1:]  # 去掉>
                current_gene = extract_gene_id(current_header)
                current_seq = []
            else:
                current_seq.append(line)
        
        # 保存最后一个序列
        if current_gene is not None:
            gene_sequences[current_gene].append({
                'header': current_header,
                'seq': current_seq
            })
    
    # 为每个基因选择最长的转录本
    selected = {}
    for gene_id, seqs in gene_sequences.items():
        if len(seqs) == 1:
            selected[gene_id] = seqs[0]
        else:
            # 选择最长的序列
            longest = max(seqs, key=lambda x: get_sequence_length(x['seq']))
            selected[gene_id] = longest
    
    # 写入输出文件
    with open(output_file, 'w') as f:
        for gene_id in sorted(selected.keys()):
            seq_data = selected[gene_id]
            if species_name:
                # 重新格式化header
                header = f">{species_name}|{seq_data['header']}"
            else:
                header = f">{seq_data['header']}"
            f.write(header + '\n')
            for i in range(0, len(seq_data['seq']), 60):
                f.write(seq_data['seq'][i] + '\n')
    
    # 统计
    total_genes = len(gene_sequences)
    single_transcript = sum(1 for seqs in gene_sequences.values() if len(seqs) == 1)
    multi_transcript = total_genes - single_transcript
    
    return {
        'total_genes': total_genes,
        'selected_sequences': len(selected),
        'single_transcript': single_transcript,
        'multi_transcript': multi_transcript
    }

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python extract_longest_pep.py <input_fasta> <output_fasta> [species_name]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    species_name = sys.argv[3] if len(sys.argv) > 3 else None
    
    print(f"处理文件: {input_file}")
    print(f"输出文件: {output_file}")
    if species_name:
        print(f"物种名: {species_name}")
    print()
    
    stats = extract_longest_transcripts(input_file, output_file, species_name)
    
    print("=" * 50)
    print("处理完成!")
    print(f"总基因数: {stats['total_genes']}")
    print(f"选择的序列数: {stats['selected_sequences']}")
    print(f"单转录本基因: {stats['single_transcript']}")
    print(f"多转录本基因: {stats['multi_transcript']} (已选择最长)")
