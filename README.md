---
title: mtdna-192component-mutspec-chordata
author: Dmitrii Iliushchenko
date: 04.11.2022
---

## Environment

- python 3.9+

## Workflow

### 1. Download sequences from RefSeq.

Using web-site [NCBI RefSeq] (https://www.ncbi.nlm.nih.gov/refseq/) with options `((complete mitochondrion RefSeq)) AND "chordates"[porgn:__txid7711] `, download GenBank file.

### 2. Get genes from GenBank

To obtain all 13 mitochondrial genes use script `scripts/get_genes_from_genbank.py`. It takes GenBank file and parses the data to 13 separated fasta files with all chordates included.
