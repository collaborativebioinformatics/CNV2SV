# CNVLINK
# CNV2SV project, BCM SV Hackathon October 2020
# M.S.- Mon, October 12 2020

#Prerequisites
#Python >= 3 (Tested on 3.8)
#conda install intervaltree
#conda install biopython
#conda install samtools

#TODO: Align unmatched CNV sequences against dup/del sequences from assembly-assembly vcf file (distant matching), this will allow checking for non-adjacent CNVs
#TODO: Confirm results

import os
import sys
import intervaltree
import collections
import subprocess
import json

version="1.0b"

cnv_vcf_file=sys.argv[1]       #CNV calls output vcf file (calls made on GRCH38 with chm13 reads)
asm_asm_vcf_file=sys.argv[2]   #dipcall output vcf file (CHM13 vs GRCH38 assembly-assembly alignment)
ref_file=sys.argv[3]           #Reference file (GRCH38)
output_dir=sys.argv[4].strip() #Output directory to store the resulting files

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

svtype_filter=["DUP"]

#Parameters for confirming CNVs where duplications / deletions are adjacent
padding_bp=1000           #Padding around CNV calls added to starts/end position when recorded in intervaltree, larger value allows for larger distance between called CNV and asm-asm event
min_cnv_length=50         #Minimum bp for a variant called in asm_asm_vcf_file to be considered as a potential duplication/deletion event
max_length_diff_frac=0.2  #Maximum length difference when matching called CNV to asm-asm events (in terms of fraction of length of longer sequence)
max_start_pos_diff=1000   #Maximum difference in starting position for adjacent duplication events
min_score_frac=0.8        #Minimum global alignment identity between CNV and asm-asm event sequence

stats={"putative_counts":{"DUP":0,"DEL":0},   #Total number of asm-asm putative dup/del events
       "cnv_call_counts":{"DUP":0,"DEL":0},   #Total number of CNV calls
       "matched_counts":{"DUP":0,"DEL":0},    #CNVs with a corresponding asm-asm event that fulfills all criteria
       "perfect_counts":{"DUP":0,"DEL":0},    #CNVs with a corresponding asm-asm event that fulfills all criteria and exactly match in terms of start position and length
       "no_overlap_counts":{"DUP":0,"DEL":0}, #CNVs with no overlapping asm-asm event
       "multihit_counts":{"DUP":0,"DEL":0},   #CNV overlaps multiple asm-asm events (currently, these CNVs are skipped!)
       "nonsense_counts":{"DUP":0,"DEL":0},   #e.g. matching DUP with DEL
       "badalign_counts":{"DUP":0,"DEL":0}}   #Start position difference too large, length difference too large or alignment too divergent between CNV sequence and asm-asm event

def readvcf(filename):
    for line_ in open(filename,"r"):
        line=line_.strip()
        if line=="" or line[0]=="#":
            continue
        chr,pos,id,ref,alt,qual,filter,info_str,format,extra=line.split("\t")[:10]
        var={"chr":chr,"pos":int(pos),"id":id,"ref":ref,"alt":alt,"qual":qual,"filter":filter,"info_str":info_str,"format":format,"extra":extra}
        yield var

def vcfinfo2dict(info_str):
    result={}
    for part in info_str.split(";"):
        if "=" in part:
            key,value=part.split("=")
            result[key]=value
    return result

def extractseq(filename,chr,start,end,header=None):
    cmd="samtools faidx %s %s:%d-%d"%(filename,chr,start,end)
    proc=subprocess.run(cmd, shell=True, capture_output=True, check=True)
    if proc.stdout==None:
        return None
    out=proc.stdout.decode("utf8")
    #n_count=sum([1 if c=="N" or c=="n" else 0 for c in out])
    #if float(n_count)/len(out) > max_n_frac:
    #    return None

    seq="".join(out.split("\n")[1:])+"\n"
    if len(seq.strip())=="":
        return None
    if header!=None:
        return ">%s"%(id)+"\n"+seq
    else:
        return seq.strip()

#STEP1: Match CNVs called on GRCH38 using short reads (cnv_vfc_file) with variant calls from assembly-assembly alignment (CHM13 vs GRCH38)

#Load potential duplications/deletions from assembly-assembly alignment variants file
trees=collections.defaultdict(lambda: intervaltree.IntervalTree())
for var in readvcf(asm_asm_vcf_file):
    length_diff=len(var["alt"])-len(var["ref"])

    if abs(length_diff) < min_cnv_length:
        continue

    putative_svtype="DUP" if length_diff>0 else "DEL"
    stats["putative_counts"][putative_svtype]+=1

    var["putative_svtype"]=putative_svtype
    var["length"]=length_diff

    trees[var["chr"]][(var["pos"]-padding_bp):(var["pos"]+padding_bp)]=var

from Bio import pairwise2
from Bio.Seq import Seq

cnvs=list(readvcf(cnv_vcf_file))

#Confirm called CNVs (when duplicated sequence is adjacent)
for var in cnvs:
    info=vcfinfo2dict(var["info_str"])
    start=var["pos"]
    end=int(info["END"])

    var["status"]="not_included"
    var["reason"]=None

    var["start"]=start
    var["end"]=end
    var["length"]=end-start
    var["svtype"]=info["SVTYPE"]
    var["info"]=info

    if svtype_filter!=None and len(svtype_filter) and not info["SVTYPE"] in svtype_filter:
        continue

    var["status"]="no_match"
    var["reason"]="none"

    stats["cnv_call_counts"][info["SVTYPE"]]+=1

    overlap=trees[var["chr"]][start:end]

    if len(overlap)==0:
        var["reason"]="no_overlap"
        stats["no_overlap_counts"][var["svtype"]]+=1
        continue

    if len(overlap)>1:
        #TODO: Implement handling of multiple hits
        var["reason"]="multihit"
        stats["multihit_counts"][var["svtype"]]+=1
        continue

    #for i,asm_var in enumerate(list(overlap)):
    asm_var=list(overlap)[0].data

    #Matched e.g. a CNV duplication call with a deletion event
    if asm_var["putative_svtype"] != var["svtype"]:
        var["reason"]="nonsense"
        stats["nonsense_counts"][var["svtype"]]+=1
        continue

    #Length difference of reported event too long
    max_length=max(var["length"],asm_var["length"])
    length_diff=abs(var["length"]-asm_var["length"])
    if max_length>0 and length_diff/float(max_length) > max_length_diff_frac:
        stats["badalign_counts"][var["svtype"]]+=1
        var["reason"]="length_diff"
        continue

    if abs(var["start"]-asm_var["pos"]) > max_start_pos_diff:
        stats["badalign_counts"][var["svtype"]]+=1
        var["reason"]="startpos_diff"
        continue

    print("cnv[",var["chr"],var["start"],var["length"],"] vs asm[",asm_var["chr"],asm_var["pos"],asm_var["length"],"]")

    asm_seq=asm_var["alt"]
    var["asm_seq"]=asm_seq

    cnv_seq=extractseq(ref_file,var["chr"],var["start"],var["end"],header=None)
    var["cnv_seq"]=cnv_seq
    var["cnv_seq_padded"]=extractseq(ref_file,var["chr"],var["start"]-30000,var["end"]+30000,header=None)
    stats["matched_counts"][var["svtype"]]+=1

    alignment=pairwise2.align.globalxx(Seq(cnv_seq), Seq(asm_seq),one_alignment_only=True)[0]
    score_frac=alignment.score/max(len(cnv_seq),len(asm_seq))
    if score_frac < min_score_frac:
        stats["badalign_counts"][var["svtype"]]+=1
        var["reason"]="alignment_score"
        continue
    var["aln_score"]=alignment.score


    if var["start"]==asm_var["pos"] and cnv_seq==asm_seq:
        stats["perfect_counts"][var["svtype"]]+=1
        var["status"]="match_adj_perfect"
    else:
        var["status"]="match_adj"

    var["matched_asm_var"]=asm_var



    #for aln in alignments:
    #    print(aln)
        #print(pairwise2.format_alignment(*aln))

    #Implemented for duplications only for now..
    #assert(putative_svtype==var["svtype"])

h=open("%s/cnvlink_adjacent.fasta"%output_dir,"w")
for var in cnvs:
    if var["status"]=="adjacent_perfect" or var["status"]=="adjacent":
        h.write(">%s\n%s\n"%(var["id"],var["cnv_seq_padded"]))
h.close()

h=open("%s/cnvlink_adjacent_out.tsv"%output_dir,"w")
h.write("#generated with cnvlink %s\n"%version)
h.write("#cmd: %s\n"%" ".join(sys.argv))
h.write("#cnv_id\tstatus\tfail_reason\ttype\tchr\taln_score\tcnv_start\tcnv_length\tasm_start\tasm_length\tcnv_seq\tasm_seq\n")
for var in sorted(cnvs,key=lambda v: v["aln_score"] if "aln_score" in v else -1,reverse=True):
    if var["status"]=="not_included" or var["status"]=="no_match":
        continue
    h.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n"%(
                                                             var["id"],
                                                             var["status"],
                                                             var["reason"],
                                                             var["svtype"],
                                                             var["chr"],
                                                             var["aln_score"] if "aln_score" in var else -1,
                                                             var["start"],
                                                             var["length"],
                                                             var["matched_asm_var"]["pos"] if "matched_asm_var" in var else -1,
                                                             var["matched_asm_var"]["length"] if "matched_asm_var" in var else -1,
                                                             var["cnv_seq"] if "cnv_seq" in var else "*",
                                                             var["asm_seq"] if "asm_seq" in var else "*"))
h.close()

print(stats)
h=open("%s/cnvlink_stats.json"%output_dir,"w")
h.write(json.dumps(stats))
h.close()
