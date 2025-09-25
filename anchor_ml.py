import os,re
import numpy as np
from Bio import SeqIO
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import classification_report, accuracy_score, f1_score
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import LabelEncoder
import pandas as pd
from sklearn.metrics import matthews_corrcoef

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
    s = {}
    fs = list(os.walk(data))[0][2]
    for filename in fs:
        if filename.endswith('.fasta') or filename.endswith('.fa'):
            filepath = os.path.join(data, filename)
            seq = SeqIO.read(filepath, 'fasta')
            key = seq.id.split('.')[0]
            s[key] = str(seq.seq)
    return s

if __name__=='__main__':
    
    # #the savepath of data, can be a file or folder
    datapath='data.fasta';
    k=4;##Parameter Settings
    epsilon=13
    
    
    if os.path.isfile(datapath):
        seqs = {i.id.split('.')[0]: str(i.seq) for i in SeqIO.parse(datapath, 'fasta')}
        
    elif os.path.isdir(datapath):
        seqs=readfolder(datapath)
        
    else:
        exit('Check the path of inputting data')
    
    print('  '+'=='*len("~"*2+'mic-alignment'+"=="*1))
    print('   '+"--"*2+'  start scanner anchor  '+"--"*2)
    
    d={seqname: scan_anchor(seq, k, epsilon) for seqname,seq in seqs.items()}
    
    allanchor=[]
    for i in d.values():
        allanchor+=list(i.keys())
    anchorcode={j:i for i,j in enumerate(set(allanchor))}
    seq_vector={}

    for i,j in d.items():
        seq_vector[i]=[0]*len(anchorcode)
        for ij,ijk in j.items():
            # let sequecne change a vector. as input vector into the ML
            seq_vector[i][anchorcode[ij]]=ijk
    
    #%% make the feature vector as a feature matrix 
    seqname_code={j:i for i,j in enumerate(d.keys())}
    ns=np.zeros((len(seqs),len(anchorcode)))        
    spacevec=[0]*len(anchorcode)   
    for i,j in d.items():
        for ij,ijk in j.items():
            ns[seqname_code[i]][anchorcode[ij]]=ijk
            
    seqnames=[0]*len(seqname_code)
    for i,j in seqname_code.items():
        seqnames[j]=i
    
#The data is my labeled data(accession and type)
metadata = pd.read_excel('data.xlsx', header=None)
sample_labels = dict(zip(metadata[0], metadata[1]))
sample_labels = {key.split('.')[0]: value for key, value in sample_labels.items()}
names = seqnames
X = ns
labels_raw = [sample_labels.get(name, None) for name in names]
mask = [label is not None for label in labels_raw]
X = X[mask]
names_filtered = np.array(names)[mask]
labels = np.array(labels_raw)[mask]
le = LabelEncoder()
y = le.fit_transform(labels)
# the test set 20%
X_train_val, X_test, y_train_val, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y)
# the five fold cross validation
n_splits = 5
skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
train_accuracies = []
val_accuracies = []
val_f1_scores = []
val_mcc_scores = []
for fold, (train_index, val_index) in enumerate(skf.split(X_train_val, y_train_val)):
    print(f"\n=== {fold+1} ===")

    X_train, X_val = X_train_val[train_index], X_train_val[val_index]
    y_train, y_val = y_train_val[train_index], y_train_val[val_index]
    
    print(f"训练集大小: {len(X_train)} ({len(X_train)/len(X_train_val):.1%})")
    print(f"验证集大小: {len(X_val)} ({len(X_val)/len(X_train_val):.1%})")
    
# the model SVM 
    model = SVC(
        kernel='rbf',
        C=6,
        gamma='scale',
        class_weight='balanced',
        decision_function_shape='ovr',
        random_state=42,
        max_iter=2000
    )
    
    # train the model
    model.fit(X_train, y_train)
    
    # accuracy
    train_accuracy = model.score(X_train, y_train)
    train_accuracies.append(train_accuracy)
    print(f"accuracy: {train_accuracy:.4f}")
    
    # in the validation set
    y_val_pred = model.predict(X_val)
    val_acc = accuracy_score(y_val, y_val_pred)
    val_f1 = f1_score(y_val, y_val_pred, average='macro')
    val_mcc = matthews_corrcoef(y_val, y_val_pred)
    val_accuracies.append(val_acc)
    val_f1_scores.append(val_f1)
    val_mcc_scores.append(val_mcc)
    print(f"accuracy={val_acc:.4f}, F1={val_f1:.4f}, MCC={val_mcc:.4f}")
# the final model
print("\n=== Final Model ===")
final_model = SVC(
    kernel='rbf',
    C=6,
    gamma='scale',
    class_weight='balanced',
    decision_function_shape='ovr',
    random_state=42,
    max_iter=2000
)
final_model.fit(X_train_val, y_train_val)
# in the test set
y_test_pred = final_model.predict(X_test)
test_acc = accuracy_score(y_test, y_test_pred)
test_f1 = f1_score(y_test, y_test_pred, average='macro')
test_mcc = matthews_corrcoef(y_test, y_test_pred)
report = classification_report(y_test, y_test_pred, target_names=le.classes_)

print("\n=== results===")
print(f"train accuracy={np.mean(val_accuracies):.4f} ± {np.std(val_accuracies):.4f}")
print(f"train F1={np.mean(val_f1_scores):.4f} ± {np.std(val_f1_scores):.4f}")
print(f"train MCC={np.mean(val_mcc_scores):.4f} ± {np.std(val_mcc_scores):.4f}")
print("\n=== test results ===")
print(f"test: {len(X_test)}")
print(f"test accuracy={test_acc:.4f}, F1分数={test_f1:.4f}, MCC={test_mcc:.4f}")
print("\nreport:")
print(report)