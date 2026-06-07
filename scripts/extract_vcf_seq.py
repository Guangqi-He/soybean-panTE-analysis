import sys
import argparse

def extract_sequences(input_vcf, sv_type, min_len, max_len, output_fasta):
    count = 0
    with open(input_vcf, 'r') as vcf, open(output_fasta, 'w') as fasta:
        for line in vcf:
            # 跳过表头
            if line.startswith('#'):
                continue
            
            cols = line.strip().split('\t')
            if len(cols) < 5:
                continue

            chrom = cols[0]
            pos = cols[1]
            ref = cols[3]
            alt = cols[4]
            
            # 处理可能有多个 ALT 的情况，只取第一个
            alt_seq = alt.split(',')[0]
            
            # 计算长度
            ref_len = len(ref)
            alt_len = len(alt_seq)
            diff_len = abs(ref_len - alt_len)

            # 长度过滤
            if not (min_len <= diff_len <= max_len):
                continue

            # 类型判断逻辑
            # INS: ALT 比 REF 长
            # DEL: REF 比 ALT 长
            current_type = ""
            var_seq = ""

            if alt_len > ref_len:
                current_type = "INS"
                # 插入序列通常是 ALT 去掉第一个锚点碱基
                var_seq = alt_seq[1:]
            elif ref_len > alt_len:
                current_type = "DEL"
                # 缺失序列通常是 REF 去掉第一个锚点碱基
                var_seq = ref[1:]
            
            # 如果类型匹配，则写入文件
            if current_type == sv_type.upper():
                fasta.write(f">{chrom}_{pos}\n")
                fasta.write(f"{var_seq}\n")
                count += 1

    print(f"Successful extracted {count}  {sv_type} sequence to {output_fasta}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="extract DEL or INS sequence from VCF file")
    parser.add_argument("-i", "--input", required=True, help="VCF file path")
    parser.add_argument("-t", "--type", required=True, choices=["DEL", "INS"], help="Variation type: DEL or INS")
    parser.add_argument("--min", type=int, default=0, help="min variation length (default: 0)")
    parser.add_argument("--max", type=int, default=100000000, help="max variation length (default: unlimited)")
    parser.add_argument("-o", "--output", required=True, help="output path")

    args = parser.parse_args()
    extract_sequences(args.input, args.type, args.min, args.max, args.output)