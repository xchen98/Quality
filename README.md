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
#### Evaluation of the quality: The genome assembly, protein sequence and annotation file can be found and downloaded on the website (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/). The Transcript comes from the publication(Zhang, R., Calixto, C.P.G., Tzioutziou, N.A., James, A.B., Simpson, C.G., Gou, W., Marquez, Y., Kalyna, M., Patro, R., Eyras, E., Barta, A., Nimmo, H.G. and Brown, J.W.S. (2015) AtRTD - A comprehensive Reference Transcript Dataset resource for accurate quantification of transcript-specific expression in Arabidopsis thaliana. New Phytologist 208, 96-101, https://ics.hutton.ac.uk/atRTD/)
 	./quality.py -a test_data/GCF_000001735.4_TAIR10.1_genomic.fna \
			    -l eukaryote_odb10 \
			    -p ./test_data/GCF_000001735.4_TAIR10.1_protein.faa \
			    -t ./test_data/AtRTD2_19April2016.fa \
			    -anno ./test_data/GCF_000001735.4_TAIR10.1_genomic.gtf \
			    -o output

#### Visualization of the results
	./visualize_multiple_assembly.py [OUTPUT_PATH] [PATH_OF_ASSEMBLIES] 
