import sys
import Bio
from Bio import Entrez
from Bio.Blast import NCBIXML
import time
start_time = time.time()
"""Returns orthologous paralogs from two .xml BLAST outputs for E-value under 10^-20"""
                         

xml1 = open(sys.argv[1])
xml2 = open(sys.argv[2])
if len(sys.argv) == 3: 
    HSP = "TRUE"
if len(sys.argv) == 4: 
    HSP = sys.argv[3]


#parse xml files
blast_records1 = NCBIXML.parse(xml1)
blast_records2 = NCBIXML.parse(xml2)


E_VALUE_THRESH = 10**-20


pairs1 = []
pairs2 = []
# iterate through xml parse
for blast_record in blast_records1:
    a = blast_record.query
    # iterate through alignments
    for alignment in blast_record.alignments:
        
        # iterate through hits and find ones below E-threshold
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                # append protein hit pairs to list in form of tuple
                b = alignment.title
                pairs1.append((a,b))
                
# do the same for second xml file                
for blast_record in blast_records2:
    a = blast_record.query
    for alignment in blast_record.alignments:
        
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                b = alignment.title
                pairs2.append((b,a))
                
# determine orthologous paralogs via set intersection
set1 = set(pairs1)
set2 = set(pairs2)
intersect = set1.intersection(set2)

# if FALSE is stated, get matches for A/B or B/A
inset2_notset1 = set2 - set1
set1_or_set2 = list(set1) + list(inset2_notset1)



# write results to file
if HSP == "TRUE":
    f = open("orthologs_nm3318_1.txt", "w")                
    for i in intersect:
        f.write(f"C. elegans: {i[0]}  -->  drosophila: {i[1]}\n")
    f.close()

# if FALSE given
if HSP == "FALSE":
    f = open("orthologs_nm3318_1_FALSE.txt", "w")                
    for i in set1_or_set2:
        f.write(f"C. elegans: {i[0]}  -->  drosophila: {i[1]}\n")
    f.close()


           



