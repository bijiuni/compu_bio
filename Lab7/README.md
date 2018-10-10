# Lab 7 - Generating Data and Compression Algorithms

## Overview

In this lab we generated some random files (binary and fasta) and grabbed some real biological information online. Compression algorithms were investigated using the data.

**NOTE: All of the data files for this analysis can be found on our [server](https://bioe131.com/user/be131-09/tree/GIT/Computational-Biology/Lab7). They were not uploaded to GitHub due to their large file size.**


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

**_A:_** It seems that in terms of level of compression, Arithmetic compress is the best for all file type we used here. Although sometimes it has the same level of performance with gzip.

<br>

**_Q:_** Which algorithm is the fastest?

**_A:_** Generally, pbzip2 is the fastedst among the algorithms we used here. In some cases, it is 20 times faster than arithmetic compress.

<br>

**_Q:_** What is the difference between bzip2 and pbzip2? Do you expect one to be faster and why?

**_A:_** pbzip2 is a modified version of bzip2. The largest difference is that pbzip2 supports multi-threading. This means that 
linear speed improvements can be achieved on multi-CPU and multi-core computers. Therefore we expect pbzip2 to be faster.

<br>

**_Q:_** How does the level of compression change as the percentage of zeros increases? Why does this
happen?

**_A:_** When the percentage of zeros increases, the level of compression improves dramatically. This is due to the intrinsic entropy of the data.

<br>

**_Q:_** What is the minimum number of bits required to store a single DNA base?

**_A:_** There are four different possibilities: A, T, G, C. The minimum number of bits required is log2(4)=2.

<br>

**_Q:_** What is the minimum number of bits required to store an amino acid letter?

**_A:_** There are 20 different possibilities. The minimum number of bits required is log2(20)=4.3. It should be an integer therefore the minimum should be 5.

<br>

**_Q:_** In your tests, how many bits did gzip and bzip2 actually require to store your random DNA and
protein sequences?

**_A:_** When it comes to storing DNA, gzip takes 29.2MB and bzip2 takes 27.3MB. Let's take 30MB for an approximation. 30MB/100,000,000 = 30 * (1,048,576 * 8)/100,000,000 bits = 2.52 bits.

When it comes to storing protein, gzip takes 63.5MB and bzip2 takes 59.8MB. Let's take 60MB for an approximation. 60/MB/100,000,000 = 60 * (1,048,576 * 8)/100,000,000 bits = 5.04 bits.

<br>

**_Q:_** Are gzip and bzip2 performing well on DNA and proteins?

**_A:_** In terms of compression ratio, bzip2 has a relatively good performance while gzip is the worst among the four compression algorithms. Arithmetic compress has the best compression rate. In terms of speed, they are both faster than arithmetic compress. However, pbzip2 is much faster than the two. In summary, they have moderate performances on DNA and proteins.
<br>


## Compressing real data


We searched using entrez using the following criteria
```
handle = Entrez.esearch(db = 'nucleotide',       # search 10 results and save the names and sequences
                        term = 'gp120 and HIV',
                        sort = 'relevance',
                        idtype = 'acc',
                        retmax = 10)
```

After checking that the results are proper, we saved the information into a dictionary:
```
dict_gp120 = dict(zip(list_name, list_seq))
```

The multi_fasta file is generated using the following methods:

```
ofile = open("multi_fasta.fa", "w")

for i in range(len(list_seq)):                    # write to the multi_fasta file
    ofile.write(">" + list_name[i] + "\n" +list_seq[i] + "\n")

ofile.close()
```

The real data is compressed using the same methods and the information can be summarized using the table below:

|original file|command type|input file size|output file size|compression ratio|
|------|------|------|------|------|
|multi_fasta.fa|gzip|6.61 kB|1.25 kB|18.91%|
|multi_fasta.fa|bzip2|6.61 kB|1.33 kB|20.12%|
|multi_fasta.fa|ArithmeticCompress	|6.61 kB|2.42 kB|36.61%|

**_Q:_**  A priori, do you expect to achieve better or worse compression here than random data? Why?

**_A:_** We expected the compression to be better than random data. As there is more pattern to be found

**_Q:_**  How does the compression ratio of this file compare to random data?

**_A:_** The compression ratio was 29.2%, 27.3%, and 25% for random data. gzip makes use of Huffman. Real sequence will have different frequency for different base. This makes sense. bzip2 uses Burrows-Wheeler. Real sequences have more patterns to be found, therefore it should perform better on real data as well. However, arithmetic compress performs worse on real data.

## Estimating compression of 1000 terabytes

The best performance we achieved was 0:00.10 for 100MB using pbzip2. For 100 TB, it would require 0.1s*10^6 = 100,000 seconds = 27.8 hours. Therefore, the key limit here is time. Furthermore, the differences between compression ratio for the four algorithms are not huge.

Assume that each data can only be processed using one computer (otherwise a large number of computing power can certainly compress all the data), we use pbzip2 for all the data.

80% re-sequencing of genomes and plasmid: One computer can compress in a day (86400 seconds/0.1 second) *100MB = 86.4TB <800TB. The reduction in data in 24 hours: (86400 seconds/0.1 second) *100MB * (1-20.12%) = 69TB

10% protein sequences: use pbzip2. One computer can compress in a day (86400 seconds/0.82 second) *100MB = 10.54TB <100TB. The reduction in data in 24 hours: (86400 seconds/0.82 second) *100MB * (1-59.8%) = 4.2TB

10% random binary data: assume this is random binary with 50% zeros. One computer can compress in a day (86400 seconds/1.5 second) *100MB = 5.76TB < 100TB. The reduction in data in 24 hours: (86400 seconds/1.5 second) *100MB * 0 = 0TB. Compression doesn't work.

Therefore, the total reduction is 49TB + 4.2TB = 53.2TB.

The money saved per day: $50/TB * 53.2TB = $2660
In a year the bonus would be 365*$2600 = $970,900
