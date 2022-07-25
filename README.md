# DenovoAS_Finder
The pipeline to annotate the teosinte transcriptome.
Usage:python DenovoAS_Finder.py
Note:First, BLAST+ and mcl should be loaded, and test.txt and data.fa should be prepared. The test.txt is the alignment result of test.fa and itself by BLAST+, which inculde 12 lines "gene label pident length mismatch gapopen qstart qend sstart send evalue bitscore". "gene" line combines isoform1 and isoform2 with *, and "label" line is consist of 0 or 1, while 0 presents that two isoforms come from different genes and 1 presents the opposite. Other lines produced by "blastn -query test.fa -db test -outfmt 6 -out test.txt -evalue 1e-10". The data.fa is the fasta file contained isoforms which you want to classify.
"""
