from scipy.spatial.distance import pdist, squareform
from functools import partial
from collections import Counter
from Bio import SeqIO

import numpy as np
import pandas as pd
import math, argparse, gc

dna={'A':0,'C':1,'G':2,'T':3}

def write_megaa(path, name, dist):
    with open(path, "w") as f:
        f.write('#mega\n!TITLE;\n')
        for namei in name:
            f.write('#' + namei + '\n')
        l = dist.shape[0]
        for i in range(l):
            for j in range(i):
                f.write(str('{:.10f}'.format(dist[i][j])) + '\t')
            f.write('\n')
    f.close()


def anchor(opt, seq, name):
    import anchor111, tqdm
    tempcheck = 1
    print('[Step 1]  obtain the   sequence (k={})'.format(opt.k))
    pa = pl()
    pas = partial(anchor111.anchor_, opt.k, opt.epsilon)
    inn_f = list(tqdm.tqdm(pa.imap_unordered(pas, seq),
                           desc='  -Progress',
                           total=len(seq), ncols=100))
    print('[Step 2]  Compute distance matricx')
    names = [];
    ns = []
    for i in inn_f:
        names.append(i[0])
        ns.append(i[1])
    del inn_f
    if tempcheck:
        distA = pdist(ns, metric=opt.metric)
        distB = squareform(distA)
    else:
        import dic_distance
        distB = dic_distance.distance_dict(ns)

    print('[Step 3]  wrirte distance matrix into file.')
    if name:
        names = [name[i] for i in names]

    savefile = opt.seqs+'_'+ str(opt.k)+'_'+str(opt.epsilon)+ '.meg'
    path = os.path.join(opt.savefolder, savefile)
    write_megaa(path, names, distB)
    # del distA, distB
    gc.collect()
    # return inn_f, distB
    return  distB


def splitn(s):
    return os.path.basename(s).replace('.fasta', '.meg')


def optil(s):
    ave = 0
    for i in range(len(s)):
        ave += len(s[i][1])
    return math.log10(ave / len(s)) - 0.5


def parameters():
    parser = argparse.ArgumentParser()

    parser.add_argument('--k', default=4, type=int, help='kmer')
    parser.add_argument('--seqs', default='hcv',
                        type=str, help='the seqs file (all sequence in one fasta file)')
    parser.add_argument('--metric', default='cosine', type=str,
                        help='metric of distance')
    parser.add_argument('--seqformat', default='fasta', type=str,
                        help='the data stype')
    parser.add_argument('--epsilon', default= 12, type=int,
                        help='the space of kmer')
    parser.add_argument('--savefolder', default='./distance', type=str,
                        help='position of save distance file')

    opt = parser.parse_args()

    return opt


if __name__ == '__main__':

    import os

    opt = parameters()
    data = os.path.join(opt.seqs)

    print('\n', '=' * 8, 'START ', '=' * 8)
    print('[Step 0]  upload data')

    if os.path.isfile(data):
        s = [];
        recordname = {}

        seq = SeqIO.parse(data, opt.seqformat)
        s = [[i.id.split('.')[0], str(i.seq)] for i in seq]

        if opt.medatacsv:
            cs = pd.read_csv(opt.medatacsv, sep='\t', header=None)
            recordname = cs.set_index(cs[0])[1].to_dict()
            recordname = {i.split('.')[0]: j for i, j in recordname.items()}
        else:
            recordname = {s[i][0]: s[i][0] for i in range(len(s))}
    elif os.path.isdir(data):
        s = [];
        recordname = {}
        fs = list(os.walk(data))[0][2]
        for i in fs:
            seq = SeqIO.read(os.path.join(data, i), opt.seqformat)
            s.append([seq.id, str(seq.seq)])
            recordname[seq.id] = i.split('.')[0]
    else:
        exit('data upload error, check your input !')

    if not os.path.exists(opt.savefolder):
        os.makedirs(opt.savefolder)
    

    print(' -k is {}'.format(opt.k))
    print(' -epsilon is {}'.format(opt.epsilon))
    print(' -total {} sequences'.format(len(s)))
    print(' -distance metric is {}'.format(opt.metric))
    print(' -save distance mega file pathe: {}'.format(
        os.path.join(opt.savefolder, opt.seqs+'_'+ str(opt.k)+'_'+str(opt.epsilon)
                     + '.meg')))
    ss = anchor(opt, s, recordname)

    print('=' * 8, 'ENDING', '=' * 8, '\n')
