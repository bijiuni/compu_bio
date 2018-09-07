# Coursework for UCB BioE 231 - Lab 2

This lab analyzes and plots the data given in the sequence file seq.fa


## Command

A step by step explanation of Lab 2

### Generating	a	phylogenetic	tree	


To align multiple sequences and generate multiple	sequence	alignment	(MSA), we used the command:

```
be231-13@jitterbug ~/lab2 [5:21pm]> muscle -in seqss.fa -out seqss.aligned.fa
```

The result:
```
MUSCLE v3.8.31 by Robert C. Edgar

http://www.drive5.com/muscle
This software is donated to the public domain.
Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.

seqss 48 seqs, max length 2217, avg  length 2209
00:00:00     11 MB(1%)  Iter   1  100.00%  K-mer dist pass 1
00:00:00     11 MB(1%)  Iter   1  100.00%  K-mer dist pass 2
00:00:07     83 MB(5%)  Iter   1  100.00%  Align node
00:00:07     83 MB(5%)  Iter   1  100.00%  Root alignment
00:00:11     83 MB(5%)  Iter   2  100.00%  Refine tree
00:00:11     83 MB(5%)  Iter   2  100.00%  Root alignment
00:00:11     83 MB(5%)  Iter   2  100.00%  Root alignment
00:00:26     83 MB(5%)  Iter   3  100.00%  Refine biparts
00:00:42     83 MB(5%)  Iter   4  100.00%  Refine biparts
00:00:42     83 MB(5%)  Iter   5  100.00%  Refine biparts
00:00:42     83 MB(5%)  Iter   5  100.00%  Refine biparts
```

Comparison of input and output of the command:

[The seqs.fa file](https://github.com/TheGaga/Lab2/blob/master/seqs.fa)
[The seqs.aligned.fa file](https://github.com/TheGaga/Lab2/blob/master/seqs.aligned.fa)


After this, we utilized Morgan	Price’s	excellent	tool	FastTree	turn	multiple	alignment	into	a	
Newick-formatted	tree with the following command:

```
be231-13@jitterbug ~/lab2 [5:41pm]> FastTree -nt < seqss.aligned.fa > trees.nwk
```

The result:
```
FastTree Version 2.1.10 SSE3
Alignment: standard input
Nucleotide distances: Jukes-Cantor Joins: balanced Support: SH-like 1000
Search: Normal +NNI +SPR (2 rounds range 10) +ML-NNI opt-each=1
TopHits: 1.00*sqrtN close=default refresh=0.80
ML Model: Jukes-Cantor, CAT approximation with 20 rate categories
Initial topology in 0.07 seconds
Refining topology: 22 rounds ME-NNIs, 2 rounds ME-SPRs, 11 rounds ML-NNIs
Total branch-length 0.954 after 0.87 sec2, 1 of 46 splits
ML-NNI round 1: LogLk = -15621.453 NNIs 8 max delta 16.59 Time 1.37
Switched to using 20 rate categories (CAT approximation)20 of 20
Rate categories were divided by 0.740 so that average rate = 1.0
CAT-based log-likelihoods may not be comparable across runs
Use -gamma for approximate but comparable Gamma(20) log-likelihoods
ML-NNI round 2: LogLk = -14251.186 NNIs 2 max delta 0.00 Time 1.69
Turning off heuristics for final round of ML NNIs (converged)
ML-NNI round 3: LogLk = -14251.086 NNIs 0 max delta 0.00 Time 2.04 (final)
Optimize all lengths: LogLk = -14251.086 Time 2.15
Total time: 2.81 seconds Unique: 48/48 Bad splits: 0/45

```

The resulting newwick will be drawn in the [iPython file](https://github.com/TheGaga/Lab2/blob/master/Lab2.ipynb).

### Identifying sequences by BLAST & Calculating sequences statistics for each cluster


All the details are explained in the [iPython file](https://github.com/TheGaga/Lab2/blob/master/Lab2.ipynb).

Here we show how the tree was clustered.

```
Cluster 1
hu.31
hu.32
hu.14
Cluster 2
hu.44
hu.46
hu.43
hu.48
Cluster 3
pi.3
pi.1
pi.2
Cluster 4
rh.43
Cluster 5
rh.58
rh.57
hu.39
rh.49
rh.51
rh.61
rh.52
rh.50
rh.53
rh.64
Cluster 6
bb.1
bb.2
rh.10
hu.17
hu.6
Cluster 7
rh.2
rh.40
hu.41
rh.38
hu.66
hu.67
hu.42
hu.37
hu.40
Cluster 8
cy.2
rh.54
rh.55
rh.48
rh.62
Cluster 9
rh.35
rh.36
rh.37
Cluster 10
cy.4
cy.6
cy.3
cy.5
rh.13
```

The cluster representatives are the first specimen of each cluster in the list. Their BLAST analysis reveals that:
* The hu.31 sequence corresponds to the _Adeno-associated virus isolate hu.31 capsid protein VP1 (cap) gene, complete dds_.
* The hu.44 sequence corresponds to the _Adeno-associated virus isolate hu.44 capsid protein VP1 (cap) gene, complete cds_.
* The pi.3 sequence corresponds to the _Adeno-associated virus isolate pi.3 capsid protein VP1 (cap) gene, complete cds_.
* The rh.58 sequence corresponds to the _Adeno-associated virus isolate rh.58 capsid protein VP1 (cap) gene, complete cds_.
* The rh.58 sequence corresponds to the _Adeno-associated virus isolate rh.58 capsid protein VP1 (cap) gene, complete cds_.
* The bb.1 sequence corresponds to the _Non-human primate Adeno-associated virus isolate AAVbb.1 capsid protein (VP1) gene, complete cds_.
* The rh.2 sequence corresponds to the _Non-human primate Adeno-associated virus isolate AAVrh.2 capsid protein (VP1) gene, complete cds_.
* The cy.2 sequence corresponds to the _Non-human primate Adeno-associated virus isolate AAVcy.2 capsid protein (VP1) gene, complete cds_.
* The rh.35 sequence corresponds to the _Non-human primate Adeno-associated virus isolate AAVrh.35 capsid protein (VP1) gene, complete cds_.
* The cy.4 sequence corresponds to the _Non-human primate Adeno-associated virus isolate AAVcy.4 capsid protein (VP1) gene, complete cds_.




## Built With

* [Matplotlib](https://matplotlib.org/) - Python 2D plotting library
* [BioPhython](https://biopython.org/) - a set of freely available tools for biological computation 


## Authors

* **Thomas Galeon** - *Initial work* - [TheGaga](https://github.com/TheGaga)
* **Zach Lyu** - *Initial work* - [bijiuni](https://github.com/bijiuni)


