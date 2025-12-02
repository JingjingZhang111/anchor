import os,re
import numpy as np
from Bio import SeqIO
import csv

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

def DNA_complement2(kmer):
    trantab = str.maketrans('ACGTacgt', 'TGCAtgca')  
    string = kmer.translate(trantab)  
    return [kmer,kmer[::-1],string[::-1]]


def bw(v,L):
    v=np.array(v)
    b=sum(list(map(lambda i: -i*np.log2(i),v[1:]/sum(v[1:]))))
    w=1-v[0]/L
    w=-w*np.log2(w)
    return b*w

def scan_anchor(sequence,k,epsilon):
    anchorset={};L=len(sequence)
    for i in range(L-k+1):
        if i+k+epsilon<len(sequence):
            kmer=sequence[i:i+k]
            if kmer in anchorset.keys():
                anchorset[kmer].append(i-anchorset[kmer][0])
                anchorset[kmer][0]=i
                continue
            if not bool(re.search(r'[^A,C,G,T]', kmer)):
                    #continue
            # anchor_test
                for j in DNA_complement2(kmer):
                    if j in sequence[i+k:i+k+epsilon]: 
                        if kmer not in anchorset.keys():
                            anchorset[kmer]=[i]
                        break
    anchorset={i:bw(j,L) for i,j in anchorset.items()}
    
    return anchorset
                    
def readfolder(data):
    s={};fs=list(os.walk(data))[0][2]
    for i in fs:
        seq=SeqIO.read(os.path.join(data,i),'fasta')
        s[seq.id]= str(seq.seq)
    return s

def opt_parameter():
    parser = argparse.ArgumentParser()

    parser.add_argument('--seqs', default='example.fasta',
                        type=str, help='the savepath of data, can be a file or folder')

    parser.add_argument('--k', default=4,
                        type=int, help='size of kmer')

    parser.add_argument('--epsilon', default=13,
                        type=int, help='The elastic area for reverse searchuence stype')

    parser.add_argument('--savename', default='results.meg',
                        type=str, help='output filename')

    opt = parser.parse_args()

    return opt


if __name__ == '__main__':
    opt = opt_parameter()

    if os.path.isfile(opt.seqs):
        seqs = {i.id: str(i.seq) for i in SeqIO.parse(opt.seqs, 'fasta')}

    elif os.path.isdir(opt.seqs):
        seqs = readfolder(opt.seqs)

    else:
        exit('Check the path of inputting data')

    print('  ' + '==' * len("~" * 2 + 'mic-alignment' + "==" * 1))
    print('   ' + "--" * 2 + '  start scanner anchor  ' + "--" * 2)

    d = {seqname: scan_anchor(seq, opt.k, opt.epsilon) for seqname, seq in seqs.items()}

    allanchor = []
    for i in d.values():
        allanchor += list(i.keys())
    anchorcode = {j: i for i, j in enumerate(set(allanchor))}
    seq_vector = {}

    for i, j in d.items():
        seq_vector[i] = [0] * len(anchorcode)
        for ij, ijk in j.items():
            seq_vector[i][anchorcode[ij]] = ijk
    seqname_code = {j: i for i, j in enumerate(d.keys())}
    ns = np.zeros((len(seqs), len(anchorcode)))
    spacevec = [0] * len(anchorcode)
    for i, j in d.items():
        for ij, ijk in j.items():
            ns[seqname_code[i]][anchorcode[ij]] = ijk

    seqnames = [0] * len(seqname_code)
    for i, j in seqname_code.items():
        seqnames[j] = i

    from scipy.spatial.distance import squareform, pdist

    distA = pdist(ns, metric='cosine')
    distB = squareform(distA)

    write_megaa(opt.savename, seqnames, distB)
