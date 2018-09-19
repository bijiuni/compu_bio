# Lab4 - SQL Database

The aim of this lab is to investigate and document the relationships between several genes, enzymes and pathways.

We chose three pathways: Glycolysis, TCA cycle, and Pentose phosphate pathway.
Four enzymes for each pathway:
Three organism: homo sapiens, drosophilia melanogaster (fruit fly), and escherichia coli K-12 MG1655


## Database design

### Genes Table
In this table we save the information of 4*3*3=36 genes.
```
Fields: id, name, description, organism, sequence
```
### Pathways table
In this table we save the information of 3 pathways.
```
Fields: id, name, description
```
### Enzymes table
In this table we save the information of 4*3*3=36 enzymes.
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

All process except the description of pathways is automated in this lab. If you want to investigate other enzymes, you can simply
edit the enzyme list at the top of the ipython file.

Entrez is the database mainly used in this lab and the most relevant result is returned in each case.

## Built With

* [sqlite3](https://www.sqlite.org/index.html) - in-process library that implements a self-contained, serverless, zero-configuration, transactional SQL database engine.
* [pandas](https://pandas.pydata.org/) - high-performance, easy-to-use data structures and data analysis tools
* [BioPhython](https://biopython.org/) - a set of freely available tools for biological computation 
* [entrez](https://www.ncbi.nlm.nih.gov/Class/MLACourse/Original8Hour/Entrez/) -  federated search engine, or web portal that allows users to search many discrete health sciences information

## Authors

* **Thomas Galeon**  - [TheGaga](https://github.com/TheGaga)
* **Zach Lyu** - [bijiuni](https://github.com/bijiuni)
