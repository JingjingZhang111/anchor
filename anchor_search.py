from itertools import product
import numpy as np
import re

def rc_seq(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    anchor_seq = ''.join(complement[base] for base in reversed(seq))
    return anchor_seq

def re_seq(seq):
    return seq[::-1]

def generate_kmer_dict(k):
    bases = 'ACGT'
    kmer_dict = {}
    for seq in product(bases, repeat=k):
        kmer = ''.join(seq)
        kmer_dict[kmer] = 0
    return kmer_dict

def kmer1(a, c, l):
    b = {}
    for k, v in a.items():
        if k[:c] not in b.keys():
            b[k[:c]] = v
        else:
            b[k[:c]] = b[k[:c]] + v

    b = {ke: list(np.sort(list(set(va)))) for ke, va in b.items()}
    gaps = {}
    for ke, va in b.items():
        if len(va) > 2:
            gap = [va[i + 1] - va[i] for i in range(len(va) - 1)]
            sd = sum(gap)
            gap = np.sum([(-i / sd) * np.log(i / sd) for i in gap])*(((l - va[-1]) / l * np.log((l - va[-1]) / l)))
        else:
            gap = 0
        gaps[ke] = gap
    return gaps

def anchor_(k, epsilon, seqs):
    name = seqs[0]
    sequence = seqs[1]
    l = len(sequence)
    anchor = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if bool(re.search(r'[^A,C,G,T]', kmer)):
            continue
        rc_kmer = rc_seq(kmer)
        re_kmer = re_seq(kmer)
        for n in range(epsilon):
            if min(i + k + n, i + k * 2 + n) > l:
                continue
            sub_seq = sequence[i + k + n:i + k * 2 + n]
            if sub_seq == rc_kmer:
                key = kmer + rc_kmer
            elif sub_seq == re_kmer:
                key = kmer + re_kmer
            elif sub_seq == kmer:
                key = kmer + kmer
            else:
                continue
            if key not in anchor:
                anchor[key] = []
            anchor[key].append(i)
    b = kmer1(anchor, k, l)
    base_dict = generate_kmer_dict(k)
    for key, value in b.items():
        if key in base_dict:
            base_dict[key] = value
    vector = [base_dict[key] for key in base_dict]
    return [name, vector]
