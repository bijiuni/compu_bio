# Lab4 - SQL Database

The aim of this lab is to investigate and document the relationships between several genes, enzymes and pathways. The data collection process is fully automated using Entrez.

3 pathways were chosen: Glycolysis, TCA cycle, and Pentose phosphate pathway.

4 enzymes were chosen for each pathway using their Enzyme Commission number (note that we chose one repetitive enzyme 4.1.1.49 in Glycolysis and TCA cycle):

| Pathway           | Enzyme 1 | Enzyme 2 | Enzyme 3 | Enzyme 4 |
|-------------------|----------|----------|----------|----------|
| Glycolysis        | 4.1.1.49  | 5.3.1.9 | 5.4.2.2 | 1.8.1.4 |
| TCA cycle         | 4.1.1.31|6.2.1.4|1.1.5.4|4.1.1.49|
| Pentose Phosphate | 3.1.3.11|2.7.1.15|3.1.1.31|5.4.2.7|

3 organisms were investigated: homo sapiens, drosophilia melanogaster (fruit fly), and escherichia coli K-12 MG1655


## Database design

### Genes Table
In this table we save the information of 4* 3 * 3 =36 genes.
```
Fields: id, name, description, organism, sequence
```
### Pathways table
In this table we save the information of 3 pathways. The description for each of these pathways was found on the [KEGG website](https://www.genome.jp/kegg/pathway.html)
```
Fields: id, name, description
```
### Enzymes table
In this table we save the information of 11 enzymes (4 * 3 - 1 repetitive). Note that we chose EC number to be the primary key.
```
Fields: id, name, description, EC
```



We designed two associative tables.

### Pathwayenzyme table
The pathway_order is the order of occurance of the enzyme in the specific pathway.
```
Fields: id, pathway, enzyme, pathway_order
```

### Enzymegene table
This table saves the relationship between the enzyems and the genes
```
Fields: id, ec, organism, gene
```

## Implementation

All processes except for the description of pathways are fully automated in this lab. To investigate other enzymes, simply edit the enzyme list at the top of the ipython file.

Entrez is the database mainly used in this lab and the most relevant result is returned in each case.

If this document does not compile with your default data rate limit, open jupyter notebook using the following command:
```
jupyter notebook --NotebookApp.iopub_data_rate_limit=10000000000
```

## Built With

* [sqlite3](https://www.sqlite.org/index.html) - in-process library that implements a self-contained, serverless, zero-configuration, transactional SQL database engine.
* [pandas](https://pandas.pydata.org/) - high-performance, easy-to-use data structures and data analysis tools
* [BioPhython](https://biopython.org/) - a set of freely available tools for biological computation 
* [entrez](https://www.ncbi.nlm.nih.gov/Class/MLACourse/Original8Hour/Entrez/) -  federated search engine, or web portal that allows users to search many discrete health sciences information

## Authors

* **Thomas Galeon**  - [TheGaga](https://github.com/TheGaga)
* **Zach Lyu** - [bijiuni](https://github.com/bijiuni)
