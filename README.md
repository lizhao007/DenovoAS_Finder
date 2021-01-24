# DenovoAS_Finder
The pipeline of annotating the teosinte transcriptome.
Usage:python DenovoAS_Finder.py
Note:First, BLAST+ and mcl should be loaded, and test.txt and data.fa should be prepared. The test.txt is the result of test.fa alignment by BLAST+, which inculde 12 lines "gene label pident length mismatch gapopen qstart qend sstart send evalue bitscore". "gene" line combines isoform1 and isoform2 with *, and "label" line is consist of 0 or 1, while 0 presents do not come from the same gene and 1 presents yes. Other lines produced by "blastn -query test.fa -db test -outfmt 6 -out test.txt -evalue 1e-10". The data.fa is the isoform file which you want to classify.
"""
