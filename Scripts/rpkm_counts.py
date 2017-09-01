# Get RPKM 

# RPKM for a given gene = 
# read counts for gene from your HTSeq analysis/
# size of gene in base pairs from GFF file, 
# for the given gene (mRNA length, incl. 5' and 3' UTRs) /
# 1000 (as we are converting to kb) / 
# total reads mapped in experiment (from XML file).

#modules=======================================================================#
import pandas as pd 
import pickle
import HTSeq
import collections
import json

#files=========================================================================#

#functions=====================================================================#

def CleanUpGtf(gtf_file):
	
	s = gtf_file[8].apply(lambda x: x.split("; "))
	gtf_file["transcript_id"] = s.apply(lambda x: x[0])
	gtf_file["gene_id"] = s.apply(lambda x: x[1])

	return gtf_file


def GTFReader(gtf):

	gtf_file = HTSeq.GFF_Reader(gtf)
	exons = HTSeq.GenomicArrayOfSets( "auto", stranded=True )

	for feature in gtf_file:
	    if feature.type == "exon":
	       exons[ feature.iv ] += feature.attr["gene_id"]
	return exons

def GetFeatures(SRR):

	counts = collections.Counter( )

	almnt_file = HTSeq.BAM_Reader(SRR)
	for almnt in almnt_file:
	   if not almnt.aligned:
	      counts[ "_unmapped" ] += 1
	      continue
	   gene_ids = set()
	   for iv, val in exons[ almnt.iv ].steps():
	      gene_ids |= val
	   if len(gene_ids) == 1:
	      gene_id = list(gene_ids)[0]
	      counts[ gene_id ] += 1
	   elif len(gene_ids) == 0:
	      counts[ "_no_feature" ] += 1
	   else:	
	      counts[ "_ambiguous" ] += 1

	# for gene_id in counts:
	#    print(gene_id, counts[ gene_id ])
	return counts

def CountsToJson(counter, file_name): 
	counterJson = json.dumps(counter)
	with open(file_name, "w") as outfile:
		json.dump(counterJson, outfile)


#constants=====================================================================#

if __name__ == "__main__":

	gtf = "/home/thuang/files/araport.gtf"
	exons = GTFReader(gtf)
	pickle_out = open("/home/thuang/files/gtf_exons.pickle", "wb")
	pickle.dump(exons, pickle_out)
	pickle_out.close()

	SRR3581889 = "/DATA/Klepikova/SRR3581889/accepted_hits.bam"
	SRR3581889 = GetFeatures(SRR3581889)
	CountsToJson(SRR3581889, "/home/thuang/files/SRR3581889_counts.json")
	# SRR3581336 = "/DATA/Klepikova/SRR3581336/accepted_hits.bam"
	# SRR3581336 = GetFeatures(SRR3581336)

	# file_name = "/home/thuang/files/SRR3581336.txt"
	# json_counter = CountsToJson(SRR3581336, file_name)

	# gtf = pd.read_csv("/home/thuang/files/araport.gtf", sep = "\t", header = None)
	# cleaned_gtf = CleanUpGtf(gtf)

	# pickle_in = open("/home/thuang/files/counts.pickle", "rb")
	# counts = pickle.load(pickle_in)
