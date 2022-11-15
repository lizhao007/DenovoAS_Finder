# Background
The pipeline is developed for classifying the transcriptome with no reference genome.

# Installation
```
git clone git@github.com:lizhao007/DenovoAS_Finder.git
```
# Requirment
BLAST+ and mcl should be installed. Other require python packages should be installed as reminder.

# Input
Test.txt and data.fa should be prepared. The test.txt is the alignment result of test.fa and itself by BLAST+, which inculde 12 lines "gene label pident length mismatch gapopen qstart qend sstart send evalue bitscore". "gene" line combines isoform1 and isoform2 with *, and "label" line is consist of 0 or 1, while 0 presents that two isoforms come from different genes and 1 presents the opposite. Other lines produced by "blastn -query test.fa -db test -outfmt 6 -out test.txt -evalue 1e-10". The data.fa is the fasta file contained isoforms which you want to classify.

# Usage
```
python DenovoAS_Finder.py
```

# Cite
Li, Z., Han, L., Luo, Z., and Li, L. Single-molecule long-read sequencing reveals extensive genomic and transcriptomic variation between maize and its wild relative teosinte (Zea mays ssp. parviglumis). Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13454
