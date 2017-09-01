# RPKM for a given gene = 
# read counts for gene from your HTSeq analysis / 
# size of gene in base pairs from GFF file, for the given gene (mRNA length, incl. 5' and 3' UTRs) / 
# 1000 (as we are converting to kb) / 
# total reads mapped in experiment (from XML file).

require(jsonlite); require(dplyr); require(tidyr); require(data.table)
require(xml2)

gtf <- readRDS("/home/thuang/files/cleaned_gtf.txt")
counts <- "/home/thuang/files/bamdata_Developmental_transcriptome.xml"
counts.json <- "/home/thuang/files/SRR3581889_counts.json"

GetTotalReadCounts <- function(counts) {
	xml <- read_xml(counts)
	cc <- xml_find_all(xml, ".//bam_file")
	SRR <- xml_attr(cc, "record_number")
	total.reads <- xml_attr(cc, "total_reads_mapped")

	reads <- as.list(total.reads)
	names(reads) <- SRR

	return(reads)
}

ExtractLongest <- function(gtf) {
	gtf %>% 
		group_by(transcript_id) %>%
		summarize(start = min(V4), end = max(V5)) -> transcripts

	transcripts %>% 
		mutate(gene_id = gsub("\\.[0123456789]+", "", transcript_id), length = end - start) %>%
		group_by(gene_id) %>%
		summarize(longest_isoform = max(length)) 
}

GetCounts <- function(counts.json) {
	jCounts = fromJSON(txt=counts.json)
}

#--------------------------------------------------------------------------------------------------------------------------------------------
reads <- GetTotalReadCounts(counts)
gene_lengths <- ExtractLongest(gtf)
jCounts <- GetCounts