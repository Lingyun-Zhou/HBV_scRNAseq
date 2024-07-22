#!/usr/bin/env python
import os
import sys
import subprocess
import re

#version
version = '0.1a'
date = '2024-06'

###### Usage information
usage = '''

    Version: %s
    Author: Fang Peng (fun_peng@foxmail.com)
    Update: %s

    Usage: %s --bam <file> --read1 <read_R1.fastq.gz> --bed <CDS/hbv_cds.bed> --out <align.txt> --sample <string>

''' % (version, date, os.path.basename(sys.argv[0]))

def hbv_bed(file):
    anno = {}
    with open(file, 'r') as f:
        full = f.readlines()
        for i in full:
            if i.strip() == '':
                continue
            infos = i.strip().split('\t')
            if infos[9] in anno:
                anno[infos[9]].append(infos)
            else:
                anno[infos[9]] = [infos]
    return(anno)

def get_Barcode_UMI(fastq_file, reads):
    if fastq_file.endswith('.gz'):
        import gzip
        fq = gzip.open(fastq_file, 'r')
    else:
        fq = open(fastq_file, 'r')
    fastq  = fq.readline().decode().strip()
    res ={}
    while fastq:
        seq = fq.readline().decode().strip()
        plus = fq.readline()
        quality = fq.readline()
        if fastq.split(' ')[0][1:] in reads:
            res[fastq.split(' ')[0][1:]] = [seq[:16], seq[16:28]]
        fastq = fq.readline().decode().strip()
    fq.close()
    return(res)

def find_near_barcode(barcode, barcode_list, max_distance = 3):
    near_barcode = {}
    for i in barcode_list:
        mis = 0
        for j in range(len(barcode)):
            if barcode[j] != i[j]:
                mis += 1
        if mis <= max_distance:
            if str(mis) in near_barcode:
                near_barcode[str(mis)].append(i)
            else:
                near_barcode[str(mis)] = [i]
    if len(near_barcode) > 0:
        near_bar = near_barcode[str(sorted([int(i)  for i in near_barcode])[0])][0]
        return(near_bar)
    else:
        return('')

def get_celltype(celltype_file, barcodes, sample):
    celltypes = {}
    with open(celltype_file, 'r') as f:
        d = f.readlines()
        for i in d:
            sample_name, barcode, celltype = i.strip().split('\t')
            if sample_name == sample:
                celltypes[barcode] = celltype
    res = {}
    for i in barcodes:
        if i in celltypes:
            res[i] = [i, celltypes[i]]
        else:
            barcode = find_near_barcode(i, celltypes.keys())
            if barcode != '':
                res[i] = [barcode, celltypes[barcode]]
            else:
                res[i] = ['', '']
    return(res)

def fix_loc(sample_name, bam, hbv_bed_file, R1_fastq, celltype_file, output):
    beds = hbv_bed(hbv_bed_file)
    step_samtools = subprocess.Popen(['samtools', 'view', '-F', '4', bam], stdout=subprocess.PIPE)
    # whole SAM line
    bamData=step_samtools.communicate()[0].decode().rstrip('\n').split('\n')
    step_samtools.stdout.close()

    final_reads = []

    for read in bamData:
        if read.strip() == '':
            continue
        recodes = read.split("\t")
        name, ref_seq, pos, cigar, seq = recodes[0], recodes[2], recodes[3], recodes[5], recodes[9]
        if cigar == "*":
            continue
        md = [i for i in recodes if 'MD:Z:' in i]
        if len(md) == 0:
            md_tag = '150'
        else:
            md_tag = md[0].replace('MD:Z:', '') 
        pattern = re.compile(r'\d+[SMID]')
        seq_pos = pattern.findall(cigar)
        left_clip_read = ''
        right_clip_read = ''
        if seq_pos[0].endswith('S'):
            softclip = seq_pos[0][:-1]
            left_clip_read = seq[:int(softclip)]
            seq = seq[int(softclip):]
        if seq_pos[-1].endswith('S'):
            softclip = seq_pos[-1][:-1]
            right_clip_read = seq[-int(softclip):]
            seq = seq[:-int(softclip):]
        if ref_seq in beds:
            start = int(beds[ref_seq][0][1]) + int(pos) - 1
            end = start + len(seq) -1
            if start <= int(beds[ref_seq][0][2]) and end > int(beds[ref_seq][0][2]) and len(beds[ref_seq]) > 1:
                break_num = end - int(beds[ref_seq][0][2])
                final_reads.append([str(i) for i in [name, beds[ref_seq][0][0], start, beds[ref_seq][0][2], seq[:-break_num], left_clip_read, '']])
                final_reads.append([str(i) for i in [name, beds[ref_seq][1][0], beds[ref_seq][1][1], int(beds[ref_seq][1][1]) + break_num - 1, seq[-break_num:], '', right_clip_read]])
            elif start > int(beds[ref_seq][0][2]):
                break_start = start - int(beds[ref_seq][0][2]) 
                break_end = end - int(beds[ref_seq][0][2]) 
                final_reads.append([str(i) for i in [name, beds[ref_seq][1][0] , int(beds[ref_seq][1][1]) + break_start - 1, int(beds[ref_seq][1][1]) + break_end - 1, seq, left_clip_read, right_clip_read]])
            else:
                final_reads.append([str(i) for i in [name, beds[ref_seq][0][0], start, end, seq, left_clip_read, right_clip_read]])

    # barcode and umi
    barcode_umi = get_Barcode_UMI(fastq_file=R1_fastq, reads={i[0] for i in final_reads})
    for i in range(len(final_reads)):
        if final_reads[i][0] in barcode_umi:
            final_reads[i] += barcode_umi[final_reads[i][0]]
        else:
            final_reads[i] += ['','']
    
    # celltype
    if celltype_file and os.path.isfile(celltype_file):
        celltypes = get_celltype(barcodes=[i[-2] for i in final_reads], celltype_file = celltype_file, sample = sample_name)
    else:
        celltypes ={}
    for i in range(len(final_reads)):
        if final_reads[i][-2] in celltypes:
            final_reads[i] += celltypes[final_reads[i][-2]]
        else:
            final_reads[i] += ['','']

    with open(output, 'w') as w:
        for i in final_reads:
            w.write('\t'.join([sample_name] + i) + '\n')

def main():
    import argparse
    ArgParser = argparse.ArgumentParser(usage = usage, formatter_class = argparse.RawTextHelpFormatter)
    ArgParser.add_argument("--version", action="version", version=version)
    ArgParser.add_argument('--bam', action="store", type=str, default=None, required=True, help="input BAM file. [%(default)s]")
    ArgParser.add_argument("--out", action="store", type=str, default=None, required=True, help="output tab file. [%(default)s]")
    ArgParser.add_argument("--read1", action="store", type=str, default=None, required=True, help="Read 1 Fastq file. [%(default)s]")
    ArgParser.add_argument("--bed", action="store", type=str, default=None, required=True, help="Input HBV bed file. [%(default)s]")
    ArgParser.add_argument("--sample", action="store", type=str, default=None, required=True, help="Sample name in celltype.  [%(default)s]")
    ArgParser.add_argument("--celltype_file", action="store", type=str, default=None, required=False, help="celltype file.  [%(default)s]")

    if len(sys.argv) == 1:
        ArgParser.parse_args('--help'.split(' '))
    args = ArgParser.parse_args()

    for file in [args.bam, args.bed, args.read1]:
        if not os.path.isfile(file):
            print("Error:\n   Can NOT find file: %s" % (file))
            sys.exit(1)
    if not os.access(os.path.dirname(os.path.abspath(args.out)), os.W_OK):
        print("Error:\n  No write permission in: %s" % (args.out))
        sys.exit(1)
    fix_loc(sample_name = args.sample, bam = args.bam, hbv_bed_file = args.bed, R1_fastq = args.read1, celltype_file = args.celltype_file, output = args.out)

if __name__ == "__main__":
    main()
