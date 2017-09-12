# RPKM for a given gene = 
# read counts for gene from your HTSeq analysis / 
# size of gene in base pairs from GFF file, for the given gene (mRNA length, incl. 5' and 3' UTRs) / 
# 1000 (as we are converting to kb) / 
# total reads mapped in experiment (from XML file).

require(jsonlite); require(dplyr); require(tidyr); require(data.table)
require(xml2)

gtf <- readRDS("~/Documents/ara/Araport11/Files/cleaned_gtf.rds")
counts <- "~/Documents/ara/Araport11/Files/bamdata_Developmental_transcriptome.xml"
counts.json <- "~/Documents/ara/Araport11/Files/SRR3581889_counts.json"

location <- "~/Documents/ara/Araport11/Files/"

GetTotalReadCounts <- function(counts) {
	xml <- read_xml(counts)
	cc <- xml_find_all(xml, ".//bam_file")
	SRR <- xml_attr(cc, "record_number")
	total.reads <- xml_attr(cc, "total_reads_mapped")

	reads <- as.list(as.integer(total.reads))
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
		summarize(longest_gene_isoform = max(length)) 
}

GetCounts <- function(counts.json) {
	jCounts <- fromJSON(txt = counts.json)
	names <- unlist(regmatches(jCounts, gregexpr("AT[0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ]+", jCounts)))
	feat <- gsub("_no_feature\\\":", "", jCounts)
	feat_ed <- gsub("_ambiguous\\\":", "", feat)
	features <- unlist(regmatches(feat_ed, gregexpr(": [0123456789]+", feat_ed))) # have _no_feature and _ambiguous
	edited_features <- substr(features, 3, nchar(features)) 
	data.frame(names, edited_features = as.integer(edited_features), stringsAsFactors = F)
}

WriteCountsandGeneLength <- function(jsonCounts, gene.lengths, file.name = "~/Documents/ara/Araport11/Files/raw_counts_and_lengths.txt") {
	doc <- left_join(jsonCounts, gene.lengths, by = c("names" = "gene_id"))
	write.table(doc, file.name, sep = "\t", quote = FALSE, row.name = FALSE)
}

GetRPKM <- function(jsonCounts, gene.lengths, reads, counts.json) {
	g.lengths <- left_join(jsonCounts, gene.lengths, by = c("names" = "gene_id"))
	bam.file <- unlist(regmatches(counts.json, gregexpr("SRR[0123456789]+", counts.json)))
	total.reads <- as.integer(reads[[bam.file]])

	g.lengths %>% 
	mutate(rpkm = edited_features/longest_gene_isoform/1000/total.reads) %>% 
	select(names, rpkm) %>%
	write.table(paste(location, paste(bam.file, "rpkm.txt", sep = "."), sep = ""), sep = "\t", row.names = FALSE, quote = FALSE) 
}


#--------------------------------------------------------------------------------------------------------------------------------------------
reads <- GetTotalReadCounts(counts)
gene.lengths <- ExtractLongest(gtf)
jsonCounts <- GetCounts(counts.json)

WriteCountsandGeneLength(jsonCounts, gene.lengths)

#SRR3581889.rpkm
rpkm <- GetRPKM(jsonCounts, gene.lengths, reads, counts.json)