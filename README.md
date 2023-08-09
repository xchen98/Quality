# Quality assessment of the genome assembly tool
 
This innovative tool offers an user-friendly approach to comprehensively assess genome assembly quality through a range of metrics based on BUSCO. Notably, it excels at identifying unannotated regions during the alignment of genome assemblies to transcriptomes or proteomes. In the end, every fragment of the genome assembly can be assigned a score based on the metric evaluation. This tool showcases versatility by allowing users to visualize results from multiple genome assemblies, facilitating easy and insightful quality comparisons.


### Requirements:
 * samtools
 * bamtocov
 * entrez direct
 * prefetch
 * fasterq-dump
 * makeblastdb
 * diamond
 * Rscript

### Usage:
#### Evaluation of the quality
 	./quality.py -a test_data/GCF_000001735.4_TAIR10.1_genomic.fna \
			    -l eukaryote_odb10 \
			    -p ./test_data/GCF_000001735.4_TAIR10.1_protein.faa \
			    -t ./test_data/AtRTD2_19April2016.fa \
			    -anno ./test_data/GCF_000001735.4_TAIR10.1_genomic.gtf \
			    -o output

#### Visualization of the results
	./visualize_multiple_assembly.py [OUTPUT_PATH] [PATH_OF_ASSEMBLIES] 
