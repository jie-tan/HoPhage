# construct Codon Markov Model
# Input: genbank file of Prokaryote in one folder
# Output: codon Markov Model of every input genebank
import argparse
import sys
import os
import numpy as np
np_load_old = np.load
np.load = lambda *a,**k: np_load_old(*a, allow_pickle=True, **k)
import pandas as pd
from Bio import SeqIO
import time
from numba import jit

parser = argparse.ArgumentParser(
    description='construct codon Markov model')
parser.add_argument('-i', dest='path_prokaryotes', nargs=1, required=True,
                    help='Path to the folder containing the genbank files of the prokaryotic genomes')
parser.add_argument('-o', dest='output_dir', nargs=1, default=[os.getcwd()],
                    help='Output directory. The default is the current path')
args = parser.parse_args()


# -i
gb_folder = os.path.abspath(os.path.expanduser(args.path_prokaryotes[0]))
if os.path.isdir(gb_folder):
    gb_files = os.listdir(gb_folder)
    gb_files.sort()
else:
    sys.exit('Input error: no such directory ' + gb_folder)

# -o
output_dir = os.path.abspath(os.path.expanduser(args.output_dir[0]))
if not os.path.isdir(output_dir):
    sys.exit('Output error: no such directory ' + output_dir)


d = dict({'A':0,'X':0,'D':0,'W':0,
          'C':1,'M':1,'S':1,'B':1,
          'G':2,'N':2,'K':2,'R':2,'V':2,
          'T':3,'H':3,'Y':3})

def str_to_num(seq,d):
    seq_num = np.zeros(len(seq),dtype=int)
    for i in range(len(seq)):
        seq_num[i] = d[seq[i]]
    return seq_num

@jit(nopython=True)
def calc(seq_num,P):
    lin = np.arange(0, int(len(seq_num) / 3) * 3, 3)
    codon = seq_num[lin] * 16 + \
            seq_num[lin + 1] * 4 + seq_num[lin + 2]
    a = codon[:-2]*64 + codon[1:-1]
    b = codon[2:]
    for n in range(len(a)):
        P[a[n], b[n]] = P[a[n], b[n]] + 1

    return P

##############################################################################
P_all = []

# gb_folder = 'refseq_archaea/'

for i in range(len(gb_files)):
    # t1 = time.clock()
    print(i)
    records = list(SeqIO.parse(gb_folder +'/'+ gb_files[i], "genbank"))
    # 统计每个gbff文件，gbff文件中可能存在多个record，比如某个原核生物还含有质粒的时候
    P = np.zeros((64 * 64, 64), dtype=int)
    for record in records:
        seq = record.seq
        seq_num = str_to_num(seq, d)
        seq_rev_num = 3 - seq_num

        flag = 0  #指示有没有cds跨越环形DNA,之前的代码没有考虑这一点，导致和之前matlab处理的数据有一点差异
        for feature in record.features:
            if feature.type == 'CDS':
                loc = feature.location
                if len(loc.parts) == 1:
                    # strand=-1 means complementary strand
                    if loc.strand == 1:
                        seq_numc = seq_num[loc.nofuzzy_start:loc.nofuzzy_end]
                    else:
                        seq_numc = seq_rev_num[loc.nofuzzy_start:loc.nofuzzy_end][::-1]
                else:
                    # 基因组有是环形的，如果一个CDS跨越了环形的连接处loc就会出现2个part
                    flag = 1
                    parts = loc.parts
                    if loc.strand == 1:
                        seq_numc = np.hstack((seq_num[parts[0].nofuzzy_start:parts[0].nofuzzy_end],seq_num[parts[1].nofuzzy_start:parts[1].nofuzzy_end]))
                    else:
                        # 这里要尤其注意一下，只有这样写是对的，改了很多次
                        seq_numc = np.hstack((seq_rev_num[parts[0].nofuzzy_start:parts[0].nofuzzy_end][::-1],seq_rev_num[parts[1].nofuzzy_start:parts[1].nofuzzy_end][::-1]))
                temp = len(seq_numc)%3
                if temp != 0:
                    # print(len(seq_numc))
                    if (seq_numc[-3:] == np.array([3, 0, 0])).all() or (seq_numc[-3:] == np.array([3, 0, 2])).all() or (
                            seq_numc[-3:] == np.array([3, 2, 0])).all():
                        # 如果结尾是终止密码子TAA TAG TGA
                        seq_numc = seq_numc[temp:]
                        # print(len(seq_numc))

                    temp = len(seq_numc) % 3
                    # 因为有发现同时满足开头是起始密码子，结束是终止密码子，应以终止密码子为准
                    if temp != 0:
                        if (seq_numc[:3] == np.array([0, 3, 2])).all() or (seq_numc[:3] == np.array([2, 3, 2])).all() or (seq_numc[:3] == np.array([3, 3, 2])).all():
                            # 如果开头是起始密码子ATG GTG TTG
                            seq_numc = seq_numc[:-temp]
                            # print(len(seq_numc))
                temp = len(seq_numc) % 3
                if temp == 0:
                    # 如果经过上面两步的矫正仍然不是3的倍数，则说明其开头结尾均不明确无法确定相位
                    # 这种无法确定的相位则不进行计算，以避免错误信息的添加
                    P = calc(seq_numc,P)
    P_all.append(P)
    # t2 = time.clock()
    # print("time:",(t2-t1))


for i in range(len(P_all)):
    P_all[i] = np.transpose(np.log(np.transpose(P_all[i] + 1) / ((np.sum(P_all[i],axis=1)) + 64)))
np.save(output_dir + '/Plog_all.npy',P_all)