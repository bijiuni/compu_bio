# Lab 7 - Generating Data and Compression Algorithms

## Overview

In this lab we generated some random files (binary and fasta) and grabbed some real biological information online. Compression algorithms were investigated using the data.

**NOTE: All of the data files for this analysis can be found on our [server](https://bioe131.com/user/be131-09/tree/GIT/Computational-Biology/Lab7). They were not uploaded to GitHub due to their large file size.**

<br>

## Background

Key information is recorded below:

data: 1000TB every day
hard disk cost: $50 per TB

Result: 1% reduction -> $500 saving per day

Options:
* Compression algorithm fast enough to compress 1000TB per day
* Not fast enough but the savings are maximized


## Simulating the data

### Binary data

Each file contains 50%-100% of zeros.

numpy functions are used to generate the data. The b here stands for binary.

```
percentages = [i for i in range(50,101,10)]  # a list a desired percentage

for p in percentages:   # call ramdom to generate random data, call packbits and then write them into file
    binary_data = np.packbits(np.random.choice([0, 1], size=8*1024*1024*100, replace=True, p = [p/100, 1-p/100]))
    open('Data/zeros_%sp' %p, 'wb').write(binary_data)
```


### DNA and protein data

Instead of saving it as a binary file, DNA and protein information is saved in two fasta files.

Note the probabilities are equal as defulat. The join function turns the numpy array into a string.

DNA:
```
dna = np.random.choice(['A', 'T', 'G', 'C'], size=100000000, replace=True);
open('Data/dna.fa', 'w').write(''.join(dna));
```

Protein:
```
protein = np.random.choice(list(string.ascii_uppercase), size=100000000, replace=True);
open('Data/protein.fa', 'w').write(''.join(protein));
```

<br>

## Compressing the data

We compressed the files using command like the following. Note that all processes are automated. We achieved this by using a
file list, the symbol !, and the symbol $.

```
!time gzip -k $file
```

More comments and codes please refer to the [ipython](https://github.com/bijiuni/compu_bio/blob/master/Lab7/Lab7.ipynb) file.

After the compression, information can be summarized into the table below:

|original file|command type|input file size|output file size|time elapse|
|------|------|------|------|------|
|zeros_50p|gzip|105 MB|105 MB|0:04.45|
|zeros_50p|bzip2|105 MB|105 MB|0:16.72|
|zeros_50p|pbzip2|105 MB|105 MB|0:01.50|
|zeros_50p|ArithmeticCompress|105 MB|105 MB|0:40.79|
|zeros_60p|gzip|105 MB|102 MB|0:05.63|
|zeros_60p|bzip2|105 MB|105 MB|0:18.21|
|zeros_60p|pbzip2|105 MB|105 MB|0:01.39|
|zeros_60p|ArithmeticCompress|105 MB|102 MB|0:42.58|
|zeros_70p|gzip|105 MB|93.6 MB|0:06.49|
|zeros_70p|bzip2|105 MB|99.8 MB|0:14.28|
|zeros_70p|pbzip2|105 MB|99.8 MB|0:01.17|
|zeros_70p|ArithmeticCompress|105 MB|92.4 MB|0:48.34|
|zeros_80p|gzip|105 MB|81.2 MB|0:16.72|
|zeros_80p|bzip2|105 MB|86.6 MB|0:12.04|
|zeros_80p|pbzip2|105 MB|86.7 MB|0:00.95|
|zeros_80p|ArithmeticCompress|105 MB|75.7 MB|0:35.39|
|zeros_90p|gzip|105 MB|58.7 MB|0:23.70|
|zeros_90p|bzip2|105 MB|61.2 MB|0:11.26|
|zeros_90p|pbzip2|105 MB|61.2 MB|0:00.76|
|zeros_90p|ArithmeticCompress|105 MB|49.2 MB|0:28.81|
|zeros_100p|gzip|105 MB|102 kB|0:00.87|
|zeros_100p|bzip2|105 MB|113 B|0:01.22|
|zeros_100p|pbzip2|105 MB|5.62 kB|0:00.10|
|zeros_100p|ArithmeticCompress|105 MB|1.03 kB|0:18.93|
|dna.fa|gzip|100 MB|29.2 MB|0:14.56|
|dna.fa|bzip2|100 MB|27.3 MB|0:09.49|
|dna.fa|pbzip2|100 MB|27.3 MB|0:00.69|
|dna.fa|ArithmeticCompress|100 MB|25 MB|0:21.41|
|protein.fa|gzip|100 MB|63.5 MB|0:04.67|
|protein.fa|bzip2|100 MB|59.8 MB|0:12.80|
|protein.fa|pbzip2|100 MB|59.8 MB|0:00.82|
|protein.fa|ArithmeticCompress|100 MB|58.8 MB|0:29.85|

<br>

## Questions

**_Q:_** Which algorithm achieves the best level of compression on each file type?

**_A:_**

<br>
**_Q:_** Which algorithm is the fastest?

**_A:_**

<br>
**_Q:_** What is the difference between bzip2 and pbzip2? Do you expect one to be faster and why?

**_A:_**

<br>
**_Q:_** How does the level of compression change as the percentage of zeros increases? Why does this
happen?

**_A:_**

<br>
**_Q:_** What is the minimum number of bits required to store a single DNA base?

**_A:_**

<br>
**_Q:_** What is the minimum number of bits required to store an amino acid letter?

**_A:_**

<br>
**_Q:_** In your tests, how many bits did gzip and bzip2 actually require to store your random DNA and
protein sequences?

**_A:_**

<br>
**_Q:_** Are gzip and bzip2 performing well on DNA and proteins?

**_A:_**
<br>


## Compressing real data


Extracting reads for a single chromosome from BAM file with samtools
```
be131-09@meowth:~/Lab6/Coverage_Plot_human$ samtools view -b human_vfast.sorted.bam chr22 > human_vfast.sorted.chr22.bam
be131-09@meowth:~/Lab6/Coverage_Plot_human$ samtools view -b human_vfast.sorted.bam chrX > human_vfast.sorted.chrX.bam
be131-09@meowth:~/Lab6/Coverage_Plot_human$ samtools view -b human_vfast.sorted.bam chrY > human_vfast.sorted.chrY.bam
```

Generate list of coordinates in chromosomes 22, X and Y and the number of times read aligned to that position
```
be131-09@meowth:~/Lab6/Coverage_Plot_human$ samtools depth -a human_vfast.sorted.chr22.bam > human_pileup_chr22.tab
be131-09@meowth:~/Lab6/Coverage_Plot_human$ samtools depth -a human_vfast.sorted.chrX.bam > human_pileup_chrX.tab
be131-09@meowth:~/Lab6/Coverage_Plot_human$ samtools depth -a human_vfast.sorted.chrY.bam > human_pileup_chrY.tab
```

### Chromosome 22
#### Coverage depth as a function of position
<img src="IMG/chr22_coverage_depth.png" width="300/">

```
# Coverage statistics being calculated
Maximum coverage =  1
Minimum coverage =  0
Average coverage = 0.007
```

#### Coverage distribution
<img src="IMG/chr22_coverage_depth_distribution.png" width="300/">

### Chromosome X
#### Coverage depth as a function of position
<img src="IMG/chrX_coverage_depth.png" width="300/">

```
Maximum coverage =  1
Minimum coverage =  0
Average coverage = 0.010
```

#### Coverage distribution
<img src="IMG/chrX_coverage_depth_distribution.png" width="300/">

### Chromosome Y
#### Coverage depth as a function of position
<img src="IMG/chrY_coverage_depth.png" width="300/">

```
Maximum coverage =  1
Minimum coverage =  0
Average coverage = 0.004
```

#### Coverage distribution
<img src="IMG/chrY_coverage_depth_distribution.png" width="300/">

### Comparing average coverage depth
The average depth for each chromosome was plotted in a bar chart. The coverage percentage was given with greater precision than in the previous iteration.

```
CALCULATING AVERAGE FOR chr22
Average depth of coverage in chr22 : 0.6868%

---------------------------------------------------------------------

CALCULATING AVERAGE FOR chrX
Average depth of coverage in chrX : 0.9715%

---------------------------------------------------------------------

CALCULATING AVERAGE FOR chrY
Average depth of coverage in chrY : 0.4118%
```

<img src="IMG/average_depth_chromosome.png" width="300/">


### Analysis
All three chromosomes have a low average coverage (less than 1%), with chromosome Y having the lowest. However, chromosome Y has the maximum coverage depth (2). This suggests some of the human DNA found in the provided sequence reads comes from an individual with a Y chromosome (rather than this match just occuring by chance). This would mean Jamie is a man. This is confirmed in the second part of our analysis (see Extra credit 2 - Zooming in on high coverage regions).

Interestingly, much like with the _S. oneidensis_ reference genome, the coverage of all three chromosomes is not uniform either. For instance, the first 15 million positions of chromosome 22 do not show any coverage.

<br>

## Extra credit 2 - Zooming in on high coverage regions
We further analysed the pileup files in our Jupyter notebook to identify regions with a higher than average coverage. Two regions of interest in chrY were identified: (18485047-18485210), (27170225-27170466). The sequences were viewed with samtools and identified using BLAST.

### chrY:18485047-18485210
```
be131-09@meowth:~/Lab6$ samtools view Coverage_Plot_human/human_vfast.sorted.bam "chrY:18485047-18485210"
read524540      0       chrY    18484909        1       300M    *       0       0       TACATTTTATGCAATCAACTTTTTTACGTGTGTGTATCTGTAGTTTCATTTCTGTGGTGATGAGTGAGACAGGTGTGGAGTAAATCAGTCCATTACATTCTTTTCTAGGTTACTTGAATTTGACATCTCAGTTCAGCATTGTAAACTCTTACAATGAACTCATAAAGTTAGAACAACGTTAAAAATAATTGCTATCTAAGTATCAGAGTTAGAATAAATTATTCCCAAGGTTTCCCCTCACTTTAAGTTTCCCTGATTCTTGTATTTTTTACTTAAATTGGATATACAATTATTATTTTT    GGEGGF;;GGGEGGGGGGG>GGGG@EGBFGGGGGF.EGGFGFBEGDGFC5GG0AAEGEG1;BGAG&GFC@CFC+2FGGBGBF-E4GGGGG>6C=GG7FGEBEGEGGEFG?GG-GGEGFGGBFFFGFGEGGGGEGFGDFFEGGBG2FEG'FFG?GGGGGGE:GGGEG5DGGGD>GG8&GFGG5GC0CGFG=+GGGGFGGEFFDGGGGGFFFFDEBCED?GGGGFFEE@?GFDFFGECDCBFGF?D%CCGGBFGGGFGG9AFF:6=GGFCFGB;5F.@D@DBD6ECGGGBFCGFBEGC4+'F    AS:i:-2 XS:i:-2 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:176A123   YT:Z:UU
read77943       16      chrY    18485047        1       300M    *       0       0       TTGTAAACTCTTACAATGAACTCATAAAGTTAGAACAAAGTTAAAAATAATTGCTATCTAAGTATCAGAGTTAGAATAAATTATTCCCAAGGTTTCCCCTCACTTTAAGTTTCCCTGATTCTTGTATTTTTTACTTAAATTGGATATACAATTACTATTTTTTCATTATTTAATTCATAATACATTTGGTAAAATAATTTCTTTTTAAGTAAAACATTTAATAGTGCAGTTTGGTTCGTGTTAATTATACTTCAACGAACCCCTTATGTTACTTGCCTAGTGACAGAGTATGTGGGTAAA    BGFGFB5=>@GFFFG'GFFGEAFFGCF4=>FFFG1D@<BD7@4E?F0DCGEGFFGGC7>:>CDG-GFGGAGGGGGEGGFB=FEGFFGFGCGEEFGEGFFGG?FFBCGAGGGEEGGGAGGDGDG@?GD/EGFGF-<:FDFGG46G9GFFAFFGFGGGFGGCGE&FGFGFGGGEGG>DFG(EEGGGFEGG:DFGGGFGGGFG:GGGGGFGFBEGAFCCFFGGFEGGGGGGGGGG7AEFEGFFGCGGG;GGGGDFFFEG:GGEEFFGGEGGGGDGEC<FGAGCGG6GFGBD4EGGGGCGGGFG    AS:i:-5 XS:i:-5 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:154T145   YT:Z:UU
```

Blasting both provided sequences yielded the same result with a 100% query cover: _Homo sapiens BAC clone RP11-455E3 from Y, complete sequence_ . This is a Bacterial Artificial Chromosome usually obtained by replication in E.Coli. From this we can infer that the gene of interest is RP11-455E3.

### chrY:27170225-27170466
```
be131-09@meowth:~/Lab6$ samtools view Coverage_Plot_human/human_vfast.sorted.bam "chrY:27170225-27170466"
read929972      0       chrY    27170165        1       300M    *       0       0       ATATGTATATCACAAGCACACTGTCAGAAAATGTCAAAGGTGAGTTTTATGATACCACACATGTCCTGTTTTTACGTGTGACAGTTGGCTGTAACCATGTGGGATGATGACAGTTATTTCTGTCAGCTGGGTTTGCATACAGGACTCACAATTTCACCTGTGTGCTGAGTGCTACATTGGTTCTGCTTGTATAACCCAAAGACTCTATAAAAGTATGTGTCAATGTTGTAATCTTTTGTGATGTTTGTACAAGAATGTGATCCATGATATCACACATGTCCCTACACCTAGTTATAAGAG    FGGGFGA4FGGGGF<%EGGFFGGG7GGGFFDCGG/GGGGGGBGGGGCDGFFEAFGADGF,,E@DGFGGGBGGGGAGGFGGG>GGGEAGFFGFGFFEGEGGGG?GEGGGGGFDGGG@GCGGGGGG=C@GF@FGGGGGDEGGFGGAAGFGGGG=GF;GE:GEGGFFF.EGGGFA@AGBGEGGGG<>GFEGGGBFGGGFEBGGGGG@AFGDGFGGGDGF/FFF9GGFGGEG><DGFD2GFBGEGG?GFEBGGGGAFGFE@2EEFFB'9G0GFGGG8/6&*>F>?GEBGDFFBFD<<G>FCGA3    AS:i:-7 XS:i:-7 XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:15T75C208 YT:Z:UU
read507912      16      chrY    27170225        1       300M    *       0       0       ATGTCCTGTTTTTACGTGTGACAGTTGGCTGCAACCATGTGGGATGATGACAGTTATTTCTGTCAGCTGGGTTTGCATACAGGACTCACAATTTCACCTGTGTGCTGAGTGCTACATTGGTTCTGCTTGTATAACCCAAAGACTCTATAAAAGTATGTGTCAATGTTGTAATCTTTTGTGATGTTTGTACAAGAATGTGATCCATGATATCACACATGTCCCTACACCTAGTTATAAGAGCCTAAATATTCTCTATTTGCTGAGTTCACATCTAAGAGTCGTTATCATTCCTGCGAGCCT    F@FB73FF0F6DB;F9GG<D:EDEG>09@-F5F<FA@BD:EDAFFFA9GGD<FGD*E>DF7FC,EGGGGD4GFFFGGF=AGGGFGGGFFFGFGF@3DBFGGGF;FFFFGFGFGGFG+0B8FGGFC5?BEGGGGGG?GCG>DFGFEFEGGGD7GBGFBFGCGDD08GGGEGFFBCGGGGGFBDCGGDFGFDGGGGGGFGGGGGGCFF/GGGGGGFGGFGGCAGGGEF=CDGEFGGGDCFGGCGGDFFEFFGGGFBGCDFGGGGG+3GCG>D*GGGF>GGGG>GGGCEFDE?CFGGGGGGFG    AS:i:0  XS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:300       YT:Z:UU
```
Blasting both provided sequences yielded the same result with a 100% query cover: _Homo sapiens chromosome Y palindromes P1, P2, P3 and inverted repeat IR2 (P1-P2-P3-IR2@) on chromosome Y_.

### Analysis
Two sequences with more than 150 nucleotides were identified in the provided reads. The length of these sequences suggests this is not a chance matching and confirms that fragments of chromosome Y DNA were indeed in the reads. Jamie is therefore a man.
