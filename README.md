Anchor
===
It is a Novel Alignment-Free Method Based on the anchor, The program is written in python.

Usage
----
python anchor_conduct.py dna --seqs $fasta --metric cosine --seqformat fasta --savefolder $savefolder --medatacsv $.csv(the sequence information) --epsilon $13 --k $4

Input Parameter
---
--k the kmer length
--epsilon type=int, help='the space of kmer'
--seqs  help='the seqs file (all sequence in one fasta file)
--metric default='cosine', type=str, help='metric of distance'
--seqformat default='fasta', type=str, help='the data stype'
--savefolder   help='position of save distance file'
--medatacsv     help='the sequence information'

Cite
---
If used，Please cite：[Runbin Tang + An effective method for identifying viral mutations based on anchors]
