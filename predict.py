import argparse
import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from numba import jit
import torch
import torch.nn as nn
from torch.utils.data import DataLoader,TensorDataset
import torch.nn.functional as F
from sklearn.preprocessing import OneHotEncoder
enc = OneHotEncoder(dtype='uint8')
from sklearn import preprocessing
min_max_scaler = preprocessing.MinMaxScaler()
import time

parser = argparse.ArgumentParser(
    description='HoPhage: an ab initio tool for identifying hosts of meta-genomic phage fragments')
parser.add_argument('-q', dest='query_phage', nargs=1, required=True,
                    help='A directory containing  all single query phage fragments with .fasta/.fna suffix OR a file with .fasta/.fna suffix which contains all query phage fragments')
parser.add_argument('-c', dest='query_phage_cds', nargs=1, required=True,
                    help='A directory containing all cds output files with .fasta/.fna suffix of single query phage fragment predicted by Prodigal OR a file with .fasta/.fna suffix which contains all cds output of query phage fragments')
parser.add_argument('-o', dest='output_dir', nargs=1, default=[os.getcwd()],
                    help='Output directory. The default is the current path')
parser.add_argument('-w', dest='weight_HoPhage_S', nargs=1, type=float, default=[0.5],
                    help='Weight of HoPhage-S. Default = 0.5')
parser.add_argument('-g',dest='genus_range',nargs=1, default=[None],
                    help='A file containing host genera names of interest')
parser.add_argument('--all',action='store_true',
                    help='If set, scores of all genera will be outputted')
args = parser.parse_args()

##############################################################################################################
# Load data

def load_data():
    print("Loading data...")
    # input1: phage fragments
    query_phage = os.path.abspath(os.path.expanduser(args.query_phage[0]))
    if os.path.isdir(query_phage):
        records_seq = []
        records_id = []
        files = os.listdir(query_phage)
        for i in range(len(files)):
            record_seq = SeqIO.parse(query_phage + '/' + files[i], "fasta")
            records_seq.append(record_seq)
            records_id.append(records_seq[i].id)

    elif os.path.isfile(query_phage):
        records_seq = list(SeqIO.parse(query_phage, "fasta"))
        records_id = []
        for i in range(len(records_seq)):
            records_id.append(records_seq[i].id)
    else:
        sys.exit('Input error: no such directory or file ' + query_phage + 'containing phage fragments')

    # input2: CDSs of phage fragments
    query_phage_cds = os.path.abspath(os.path.expanduser(args.query_phage_cds[0]))
    if os.path.isdir(query_phage_cds):
        records_cds = [[] for i in range(len(records_seq))]
        files = os.listdir(query_phage_cds)
        for i in range(len(files)):
            record_cds = SeqIO.parse(query_phage_cds + '/' + files[i], "fasta")
            temp = record_cds.id.split('_')
            id = "_".join(temp[:-1])  # 因为ID中也可能存在_ 所以分割之后要再把id拼接一下
            if id in records_id:
                ind = records_id.index(id)
                records_cds[ind].append(record_cds)

    elif os.path.isfile(query_phage_cds):
        records_cds0 = list(SeqIO.parse(query_phage_cds, "fasta"))
        records_cds = [[] for i in range(len(records_seq))]
        for i in range(len(records_cds0)):
            temp = records_cds0[i].id.split('_')
            id = "_".join(temp[:-1])  # 因为ID中也可能存在_ 所以分割之后要再把id拼接一下
            if id in records_id:
                ind = records_id.index(id)
                records_cds[ind].append(records_cds0[i])
    else:
        sys.exit('Input error: no such directory or file ' + query_phage + 'containing CDSs of phage fragments')

    # ouput
    output_dir = os.path.abspath(os.path.expanduser(args.output_dir[0]))
    if not os.path.isdir(output_dir):
        sys.exit('Output error: no such directory ' + output_dir)

    # weight of HoPhage-S
    weight_s = args.weight_HoPhage_S[0]
    if weight_s is not None:
        if ((weight_s>1) | (weight_s<0)):
            sys.exit('Values error: the weight of HoPhage-S should be a value between 0 and 1')

    # genus range
    genus_interest = args.genus_range[0]
    if genus_interest is not None:
        if not os.path.isfile(genus_interest):
            sys.exit('Input error: no such file ' + genus_interest + 'host genera names of interest')
        else:
            genus_range = pd.read_csv(genus_interest)
            genus_range = genus_range['Genus'].values
    else:
        genus_range = genus_interest

    return records_id,records_seq,records_cds,output_dir,weight_s,genus_range

##############################################################################################################
# data preprocessing

d = dict({'A':0,'X':0,'D':0,'W':0,
          'C':1,'M':1,'S':1,'B':1,
          'G':2,'N':2,'K':2,'R':2,'V':2,
          'T':3,'H':3,'Y':3})

def str_to_num(seq,d):
    seq_num = np.zeros(len(seq),dtype='uint8')
    for i in range(len(seq)):
        seq_num[i] = d[seq[i]]
    return seq_num

def preprocessing(records_seq,records_cds):
    print("Data preprocessing...")
    # input
    # records_seq: Original sequence data of all query phage fragments
    # records_cds: cds prediction data from Prodigal of all query phage fragments
    # output
    # SEQ_all: Convert the original sequence to a number format
    # CDS_all: Convert the sequence of CDS region to a number format
    # SEQ_01: 1 means CDS region, 0 means not
    SEQ_all = [[] for i in range(len(records_seq))]
    CDS_all = [[] for i in range(len(records_seq))]
    SEQ_01 = [[] for i in range(len(records_seq))]

    for i in range(len(records_seq)):
        seq = records_seq[i].seq
        seq_num = np.array(str_to_num(seq, d))
        seq_rev_num = 3 - seq_num
        SEQ_all[i] = np.array([seq_num,seq_rev_num[::-1]])

        CDS = []
        seq_01 = np.zeros((2, len(str(records_seq[i].seq))),dtype='uint8')
        for j in range(len(records_cds[i])):
            cds = records_cds[i][j].seq
            cds_num = np.array(str_to_num(cds, d))
            CDS.append(cds_num)

            cds_description = records_cds[i][j].description
            ind1 = list((pos for pos, val in enumerate(cds_description) if (val == '#')))
            cds_start = int(cds_description[ind1[0] + 2: ind1[1] - 1])
            cds_end = int(cds_description[ind1[1] + 2: ind1[2] - 1])
            strand = int(cds_description[ind1[2] + 2: ind1[3] - 1])
            if strand == 1:
                seq_01[0, cds_start - 1:cds_end] = 1
            else:
                seq_01[1, cds_start - 1:cds_end] = 1

        CDS_all[i] = CDS
        seq_01[1,:] = seq_01[1,:][::-1]
        SEQ_01[i] = seq_01

    return SEQ_all,CDS_all,SEQ_01

##############################################################################################################
# function for HoPhage-S
@jit(nopython=True)
def cals(seq_num,score,P_all):
    # lengths of some cds are not multiples of 3
    lin = np.arange(0, int(len(seq_num) / 3) * 3, 3)
    codon = seq_num[lin] * 16 + \
            seq_num[lin + 1] * 4 + seq_num[lin + 2]
    a = codon[0:-2] * 64 + codon[1:-1]
    b = codon[2:]
    c = len(a)
    for n in range(c):
        score = score + P_all[:,a[n],b[n]]
    return score,c

def mkscore(cds_yn,P_all,CDS):
    if cds_yn == 1:
        score = np.zeros(np.shape(P_all)[0])
        codon_num = []
        for j in range(len(CDS)):
            seq_num = CDS[j]
            score, c = cals(seq_num, score, P_all)
            codon_num.append(c)
        score = score / sum(codon_num)
    else:
        score_temp = []
        # 3 phases in forward and reverse sequence, 6 phases in total
        for m in range(2):
            seq_num = CDS[m,:]
            for j in range(3):
                score = np.zeros(np.shape(P_all)[0])
                score, c = cals(seq_num[j:len(seq_num)], score, P_all)
                score_temp.append(score / c)
        score = np.array(score_temp).max(axis=0)
    return score

##############################################################################################################
# function for HoPhage-G

class Inception(nn.Module):
    def __init__(self, in_ch, out_ch1, mid_ch13, out_ch13, mid_ch15, out_ch15, out_ch_pool_conv):
        super(Inception, self).__init__()

        self.conv1 = nn.Sequential(
            nn.Conv1d(in_ch, out_ch1, kernel_size=1, stride=1),
            nn.ReLU())
        self.conv13 = nn.Sequential(
            nn.Conv1d(in_ch, mid_ch13, kernel_size=1, stride=1),
            nn.ReLU(),
            nn.Conv1d(mid_ch13, out_ch13, kernel_size=3, stride=1, padding=1),
            nn.ReLU())

        self.conv15 = nn.Sequential(
            nn.Conv1d(in_ch, mid_ch15, kernel_size=1, stride=1),
            nn.ReLU(),
            nn.Conv1d(mid_ch15, out_ch15, kernel_size=5, stride=1, padding=2),
            nn.ReLU())

        self.pool_conv1 = nn.Sequential(
            nn.MaxPool1d(3, stride=1, padding=1),
            nn.Conv1d(in_ch, out_ch_pool_conv, kernel_size=1, stride=1),
            nn.ReLU())

    def forward(self, inputs, train=False):
        conv1_out = self.conv1(inputs)
        conv13_out = self.conv13(inputs)
        conv15_out = self.conv15(inputs)
        pool_conv_out = self.pool_conv1(inputs)
        outputs = torch.cat([conv1_out, conv13_out, conv15_out, pool_conv_out], 1)  # depth-wise concat

        return outputs

class MyNet(nn.Module):
    def __init__(self):
        super(MyNet, self).__init__()

        self.layer12_0 = nn.Sequential(
            nn.Conv1d(6, 64, 7, 2, 3),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.MaxPool1d(3, 2, 1),
            nn.Conv1d(64, 64, 1),
            nn.ReLU(),
            nn.Conv1d(64, 192, 3, 1, 1),
            nn.BatchNorm1d(192),
            nn.ReLU(),
            nn.MaxPool1d(3, 2, 1)
        )
        # in_ch,out_ch_1,mid_ch_13,out_ch_13,mid_ch_15,out_ch_15,out_ch_pool_conv
        self.layer12_1 = nn.Sequential(
            Inception(192, 64, 96, 128, 16, 32, 32),
            Inception(256, 128, 128, 256, 32, 64, 64),
            # nn.AvgPool1d(3, 2, 1)
            nn.AdaptiveAvgPool1d(1)
        )

        self.layer3_0 = nn.Sequential(
            nn.Conv1d(in_channels=64, out_channels=512, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            # nn.AvgPool1d(kernel_size=3, stride=2, padding=1),
            nn.AdaptiveAvgPool1d(1)
        )

        self.layer4_0 = nn.Sequential(
            nn.Linear(1024,512),
            nn.ReLU(),
            nn.Dropout(0.3)
        )

        self.layer5 = nn.Sequential(
            nn.Linear(512*3, 512),
            nn.ReLU(),
            nn.Dropout(0.3)
        )

        self.output = nn.Linear(512, 2)
        # self.output = nn.Linear(512, 1)

    # forward
    def forward(self, input12,input3,input4):
        out12 = self.layer12_0(input12)
        out12 = self.layer12_1(out12)

        out3 = self.layer3_0(input3)

        out4 = self.layer4_0(input4)

        # concat
        out = torch.cat([out12.view(-1,512), out3.view(-1,512), out4.view(-1,512)], 1)
        out = self.layer5(out)
        out = self.output(out)

        return out

def pred(loader, net, device):
    # change to eval mode during testing or prediction
    net.to(device)
    net.eval()
    predicted_result = []
    for i, data in enumerate(loader, 0):
        input12, input3, input4 = data
        input12 = input12.type(torch.FloatTensor)
        input3 = input3.type(torch.FloatTensor)
        input4 = input4.type(torch.FloatTensor)
        input12, input3, input4 = input12.to(device), input3.to(device), input4.to(device)
        with torch.no_grad():
            outputs = net(input12, input3, input4)
        if i==0:
            predicted_result = outputs.cpu().numpy()
        else:
            predicted_result = np.concatenate((predicted_result,outputs.cpu().numpy()),axis=0)
    temp = F.softmax(torch.from_numpy(predicted_result), dim=1).numpy()
    result = temp[:, 1]
    return result

##############################################################################################################
# Load model

def load_model(ind_dl,ind_mk,device):
    print("Loading model...")

    # HoPhage-G
    net1 = MyNet()
    net2 = MyNet()
    net3 = MyNet()

    if device.type == 'cuda':
        print("device: %s" % device)
        # load 到GPU上
        checkpoint = torch.load('model/model_100-400.pth')
        net1.load_state_dict(checkpoint['model_state_dict'])
        checkpoint = torch.load('model/model_400-800.pth')
        net2.load_state_dict(checkpoint['model_state_dict'])
        checkpoint = torch.load('model/model_800-1200.pth')
        net3.load_state_dict(checkpoint['model_state_dict'])
    else:
        # load到CPU上
        print("device: %s" % device)
        checkpoint = torch.load('model/model_100-400.pth',
                                map_location=lambda storage, loc: storage)
        net1.load_state_dict(checkpoint['model_state_dict'])
        checkpoint = torch.load('model/model_400-800.pth',
                                map_location=lambda storage, loc: storage)
        net2.load_state_dict(checkpoint['model_state_dict'])
        checkpoint = torch.load('model/model_800-1200.pth',
                                map_location=lambda storage, loc: storage)
        net3.load_state_dict(checkpoint['model_state_dict'])

    DICODON = np.load("model/DICODON.npy")[ind_dl]
    MER = np.load("model/MER.npy")[ind_dl]

    # HoPhage-S
    MK_host = np.load('model/Plog_all.npy')[ind_mk]

    return net1,net2,net3,DICODON,MER,MK_host



if __name__ == '__main__':
    # t0 = time.time()
    records_id,records_seq,records_cds,output_dir,weight_s,genus_range = load_data()
    SEQ_all, CDS_all, SEQ_01 = preprocessing(records_seq,records_cds)
    host_taxa = pd.read_csv("host_taxa_1353.csv")
    host_info_a = pd.read_csv("info_refseq_archaea.csv", delimiter='\t', keep_default_na=False)
    host_info_b = pd.read_csv("info_refseq_bacteria.csv", delimiter='\t', keep_default_na=False)
    host_info = pd.concat([host_info_a, host_info_b], ignore_index=True)
    if genus_range is not None:
        ind_dl = host_taxa[host_taxa['genus'].isin(genus_range)].index
        host_taxa = host_taxa.iloc[ind_dl].reset_index()
        ind_mk = host_info[host_info['genus'].isin(genus_range)].index
        host_info = host_info.iloc[ind_mk].reset_index()
    else:
        ind_dl = np.arange(0,len(host_taxa))
        ind_mk = np.arange(0, len(host_info))

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    # device = torch.device("cpu" if torch.cuda.is_available() else "cpu")
    net1, net2, net3, DICODON, MER, MK_host = load_model(ind_dl,ind_mk,device)

    # t1 = time.time()
    # print("time:",(t1-t0))
    ####################################################################################
    # Score by HoPhage-G
    print("Scoring by HoPhage-G...")

    # candidate host
    DICODON = torch.from_numpy(DICODON)
    MER = torch.from_numpy(MER)
    MER = MER.unsqueeze(1)
    input3 = DICODON.type(torch.FloatTensor)
    input4 = MER.type(torch.FloatTensor)

    Score_dl = []
    for i in range(len(SEQ_all)):
        # print(i)
        SEQ = SEQ_all[i]
        SEQ01 = SEQ_01[i]
        len_seq = np.shape(SEQ)[1]
        if len_seq > 1200:
            # lenght of phage fragment longer than 1200bp
            SEQ[1, :] = SEQ[1, :][::-1]
            SEQ01[1, :] = SEQ01[1, :][::-1]
            num1 = len_seq // 1200
            num2 = len_seq % 1200
            W = []
            NET = []
            FRAG = []
            FRAG01 = []
            for j in range(num1):
                W.append(1200)
                NET.append(net3)
                frag = SEQ[:, j * 1200:(j + 1) * 1200]
                frag[1, :] = frag[1, :][::-1]
                FRAG.append(frag)
                frag01 = SEQ01[:, j * 1200:(j + 1) * 1200]
                frag01[1, :] = frag01[1, :][::-1]
                FRAG01.append(frag01)
            # only the last frame which length longer than 400bp will kept
            if num2 >= 400:
                W.append(num2)
                frag = SEQ[:, num1 * 1200:len_seq]
                frag[1, :] = frag[1, :][::-1]
                frag01 = SEQ01[:, num1 * 1200:len_seq]
                frag01[1, :] = frag01[1, :][::-1]
                if num2 > 800:
                    NET.append(net3)
                    frag = np.pad(frag, ((0, 0), (0, 1200 - num2)), 'constant', constant_values=4)
                    frag01 = np.pad(frag01, ((0, 0), (0, 1200 - num2)), 'constant', constant_values=2)
                else:
                    NET.append(net2)
                    frag = np.pad(frag, ((0, 0), (0, 800 - num2)), 'constant', constant_values=4)
                    frag01 = np.pad(frag01, ((0, 0), (0, 800 - num2)), 'constant', constant_values=2)
                FRAG.append(frag)
                FRAG01.append(frag01)
            score_temp = []
            for j in range(len(FRAG)):
                enc.fit([[0], [1], [2], [3], [4]])
                onehot5 = enc.transform(np.hstack((FRAG[j][0], FRAG[j][1]))[:, np.newaxis]).toarray()
                seq = torch.from_numpy(onehot5[:, 0:4])
                enc.fit([[0], [1], [2]])
                onehot3 = enc.transform(np.hstack((FRAG01[j][0], FRAG01[j][1]))[:, np.newaxis]).toarray()
                seq01 = torch.from_numpy(onehot3[:, 0:2])

                input12 = torch.cat((seq, seq01), dim=1)
                input12 = input12.unsqueeze(0)
                input12 = input12.permute(0, 2, 1)
                input12 = input12.repeat(len(DICODON), 1, 1)
                input12 = input12.type(torch.FloatTensor)

                testdataset = TensorDataset(input12, input3, input4)
                testloader = DataLoader(testdataset, batch_size=512, shuffle=False, num_workers=0)
                # prediction
                score_dl = pred(testloader, NET[j], device)
                score_temp.append(score_dl)
            score_temp = np.array(score_temp)
            W = np.array(W).reshape(len(W), 1)
            score_dl = sum(score_temp * W) / sum(W)
            Score_dl.append(score_dl)

        else:
            # lenght of phage fragment shorter than 1200bp
            if len_seq > 800:
                len_frag = 1200
                net = net3
            else:
                if len_seq > 400:
                    len_frag = 800
                    net = net2
                else:
                    len_frag = 400
                    net = net1

            SEQ = np.pad(SEQ, ((0, 0), (0, len_frag - len_seq)), 'constant', constant_values=4)
            enc.fit([[0], [1], [2], [3], [4]])
            onehot5 = enc.transform(np.hstack((SEQ[0], SEQ[1]))[:, np.newaxis]).toarray()
            seq = torch.from_numpy(onehot5[:, 0:4])
            SEQ01 = np.pad(SEQ01, ((0, 0), (0, len_frag - len_seq)), 'constant', constant_values=2)
            enc.fit([[0], [1], [2]])
            onehot3 = enc.transform(np.hstack((SEQ01[0], SEQ01[1]))[:, np.newaxis]).toarray()
            seq01 = torch.from_numpy(onehot3[:, 0:2])

            input12 = torch.cat((seq, seq01), dim=1)
            input12 = input12.unsqueeze(0)
            input12 = input12.permute(0, 2, 1)
            input12 = input12.repeat(len(DICODON), 1, 1)
            input12 = input12.type(torch.FloatTensor)

            testdataset = TensorDataset(input12, input3, input4)
            testloader = DataLoader(testdataset, batch_size=512, shuffle=False, num_workers=0)
            # 测试
            score_dl = pred(testloader, net, device)
            Score_dl.append(score_dl)
    Score_dl = np.array(Score_dl)

    # t2 = time.time()
    # print("G time:", (t2 - t1))
    ####################################################################################
    # Score by HoPhage-S
    print("Scoring by HoPhage-S...")

    Score_mk = []
    for i in range(len(SEQ_all)):
        # print(i)
        CDS = CDS_all[i]
        SEQ = SEQ_all[i]
        cds_yn = 1
        if len(CDS) == 0:
            cds_yn = 0   # indicating that there is no cds region in this phage frgments are predited by Prodigal
            CDS = SEQ
        score = mkscore(cds_yn, MK_host, CDS)
        Score_mk.append(score)
    Score_mk = np.array(Score_mk)

    # t3 = time.time()
    # print("S time:", (t3 - t2))
    ####################################################################################
    # Integrate results from HoPhage-G and HoPhage-S
    print("Integrating results...")

    # Only top 1 genus are outputted
    output_top1 = pd.DataFrame(columns=['ID','Score-G','Score-S','Integrated Score','Host Name','Superkingdom','Phylum','Class','Order','Family','Genus'])
    for i in range(len(Score_dl)):
        score_dl = Score_dl[i]
        ind1 = np.argsort(-score_dl)
        score_dl = sorted(score_dl, reverse=True)

        ind2 = np.where(np.array(score_dl) > 0.8)[0]
        if len(ind2) > 0:
            # Maximum score of HoPhage-G is greater than 0.8
            # w = 0.125
            genus_top = host_taxa['genus'][ind1[ind2]].values
            ind_top = host_info[host_info['genus'].isin(genus_top)].index
        else:
            ind3 = np.where(np.array(score_dl) > 0.4)[0]
            if len(ind3) > 0:
                # Maximum score of HoPhage-G is between 0.4 and 0.8
                # w = 0.25
                ind4 = np.where(np.array(score_dl) > 0.25)[0]
                genus_top = host_taxa['genus'][ind1[ind4]].values
                ind_top = host_info[host_info['genus'].isin(genus_top)].index
            else:
                # Maximum score of HoPhage-G is less than 0.4
                # w = 0.5
                genus_top = host_taxa['genus'][ind1].values
                ind_top = np.arange(0,len(host_info))

        score_mk = Score_mk[i, ind_top]

        score_mk2 = np.array(sorted(score_mk, reverse=True)).reshape(-1, 1)
        # Map the score between 0 and the maximum score of A
        # score_mk3 = min_max_scaler.fit_transform(score_mk2)*np.max([score_dl[0],0.4])
        score_mk3 = min_max_scaler.fit_transform(score_mk2) * score_dl[0]
        ind5 = np.argsort(-score_mk)
        genus_sortmk, i_sortmk = np.unique(host_info['genus'][ind_top[ind5]].values, return_index=True)
        score = np.zeros((1, len(genus_top)))

        for n in range(len(genus_top)):
            score[0, n] = (1 - weight_s) * score_dl[n] + weight_s * \
                          score_mk3[i_sortmk[np.where(genus_sortmk == genus_top[n])[0][0]]][0]

        ind7 = np.where(score == np.max(score))[1][0]
        pre_genus = genus_top[ind7]
        pre_taxa = host_taxa[host_taxa['genus'] == pre_genus]

        ind8 = ind5[i_sortmk[np.where(genus_sortmk == genus_top[ind7])[0][0]]]
        ind9 = ind_top[ind8]
        info = {'ID':records_id[i], 'Score-G':score_dl[ind7], 'Score-S':score_mk[ind8], 'Integrated Score':score[0,ind7],
                'Host Name': host_info['organism_name'][ind9], 'Superkingdom': host_info['superkingdom'][ind9], 'Phylum':host_info['phylum'][ind9],
                'Class':host_info['class'][ind9], 'Order':host_info['order'][ind9],
                'Family':host_info['family'][ind9], 'Genus':host_info['genus'][ind9]}
        output_top1 = output_top1.append(info,ignore_index=True)
    # Save
    output_top1.to_csv(output_dir + '/' + "top1_prediction.csv",index=False)

    # All genera are outputted
    if args.all:
        print("All genera will be outputted")
        for i in range(len(Score_dl)):
            score_dl = Score_dl[i]
            ind1 = np.argsort(-score_dl)
            score_dl = sorted(score_dl, reverse=True)
            genus_top = host_taxa['genus'][ind1].values
            ind_top = np.arange(0, len(host_info))

            score_mk = Score_mk[i]
            score_mk2 = np.array(sorted(score_mk, reverse=True)).reshape(-1, 1)
            # Map the score between 0 and the maximum score of A
            score_mk3 = min_max_scaler.fit_transform(score_mk2) * np.max([score_dl[0],0.5])
            ind5 = np.argsort(-score_mk)
            genus_sortmk, i_sortmk = np.unique(host_info['genus'][ind_top[ind5]].values, return_index=True)
            score = np.zeros((1, len(genus_top)))

            for n in range(len(genus_top)):
                score[0, n] = (1 - weight_s) * score_dl[n] + weight_s * \
                              score_mk3[i_sortmk[np.where(genus_sortmk == genus_top[n])[0][0]]][0]

            ind7 = np.argsort(-score)
            output_temp = pd.DataFrame(columns=['ID','Score-G','Score-S','Integrated Score','Host Name','Superkingdom','Phylum','Class','Order','Family','Genus'])
            for n in range(len(ind7[0])):
                ind8 = ind5[i_sortmk[np.where(genus_sortmk == genus_top[ind7[0,n]])[0][0]]]
                ind9 = ind_top[ind8]
                info = {'ID': records_id[i], 'Score-G': score_dl[ind7[0,n]], 'Score-S': score_mk[ind8],
                        'Integrated Score': score[0, ind7[0,n]],
                        'Host Name': host_info['organism_name'][ind9], 'Superkingdom': host_info['superkingdom'][ind9],
                        'Phylum': host_info['phylum'][ind9],
                        'Class': host_info['class'][ind9], 'Order': host_info['order'][ind9],
                        'Family': host_info['family'][ind9], 'Genus': host_info['genus'][ind9]}
                output_temp = output_temp.append(info, ignore_index=True)
            output_temp.to_csv(output_dir + '/' + records_id[i] + "_prediction.csv", index=False)

    ####################################################################################
    print("Completed!")
    # t4 = time.time()
    # print("all time:", (t4 - t0))
